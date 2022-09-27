function GMMStruct_KDG = KEM_initGuess(x,h,M,MaxIter,tol,GMMStruct_0)

[N,p] = size(x);

GMMStruct_prime = GMMStruct_0;

GMMStruct = GMMStruct_prime;

omega = calc_w(x,h);

% Start Algorithm %

iteration_step_outter = 1;

crit = 10^6;

while( crit >= tol && iteration_step_outter < MaxIter) 

    % E - Step %

    gamma_prime = calc_gamma(x,GMMStruct_prime);

    for m = 1:M
      
        GMMStruct.Alpha(m) = gamma_prime(:,m)' * omega;

        GMMStruct.Means(:,m) = ( sum(gamma_prime(:,m).*x) + 1e-200 ) / ( sum(gamma_prime(:,m)) + 1e-200 );

        dist = (x-GMMStruct_prime.Means(:,m)') .* ...
            sqrt( ( gamma_prime(:,m) + 1e-200 ) / ( sum(gamma_prime(:,m)) + 1e-200 ) );

        GMMStruct.Covars(:,:,m) = dist' * dist;

    end
    
    GMMStruct.Alpha = GMMStruct.Alpha ./ sum(GMMStruct.Alpha);

    for r = 1:M
            
            curr_det = det(GMMStruct.Covars(:,:,m));
            curr_cond = cond(GMMStruct.Covars(:,:,m));
            pos_eig = all(eig(GMMStruct.Covars(:,:,m)) > 0);

            if curr_cond > 300 || abs(curr_det) > 10^3 || ~ pos_eig
               
                iteration_step_outter = MaxIter;
                
                GMMStruct = GMMStruct_prime;

            end
            
    end
    
    GMMStruct.Covars = fix_covars(GMMStruct.Covars,M);
    
    iteration_step_outter = iteration_step_outter + 1;
    
    crit = calc_loglike(x,N,GMMStruct,GMMStruct_prime);

    GMMStruct_prime = GMMStruct;

end

GMMStruct_KDG = GMMStruct;