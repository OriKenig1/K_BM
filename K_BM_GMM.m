function [GMMStruct_K,h] = K_BM_GMM(x,M,h,GMMStruct_0,MaxIter,tol)

[N,~] = size(x);

% Initialize %

GMMStruct = KEM_initGuess(x,h,M,MaxIter,tol,GMMStruct_0);

GMMStruct_prime = GMMStruct;

omega = calc_w(x,h);

% Algorithm Start %

iteration_step_outter = 1;

crit_outter = 10^6;

while( crit_outter >= tol && iteration_step_outter < MaxIter) 

    % Start of B - Step %

    u_prime = calc_u(x,GMMStruct_prime,h);

    gamma_prime = calc_gamma(x,GMMStruct_prime);

    lambda_prime = calc_lambda(gamma_prime,omega);
   
    % End of B - Step %

    % Start of M - Step %

    iteration_step_inner = 1;

    crit_inner = 10^6;

    GMMStruct_prime2 = GMMStruct_prime;

    while( crit_inner >= tol && iteration_step_inner < MaxIter) 
        
        % Calculate support parameters %

        u = calc_u(x,GMMStruct_prime2,h);

        eta = u / ( u_prime + 1e-200 );

        xi = calc_xi(x, GMMStruct, h, u_prime);

        rho = xi .* GMMStruct_prime2.Alpha';

        % Calculate main parameters %

        GMMStruct.Alpha = calc_alpha(lambda_prime, eta, xi, GMMStruct_prime2.Alpha);

        GMMStruct.Means = calc_means(x,h,GMMStruct_prime2,rho,lambda_prime,omega,gamma_prime);

        [GMMStruct.Covars,exit] = calc_covars(x,h,GMMStruct_prime2,omega,gamma_prime,rho,lambda_prime);

        % Finish inner iteration

        iteration_step_inner = iteration_step_inner + 1;

        crit_inner = calc_loglike(x,N,GMMStruct,GMMStruct_prime2);

        GMMStruct_prime2 = GMMStruct;

        if exit
            
            iteration_step_inner = MaxIter;
            iteration_step_outter = MaxIter;

        end

    end
    
    % End of M - Step %

    iteration_step_outter = iteration_step_outter + 1;

    crit_outter = calc_loglike(x,N,GMMStruct,GMMStruct_prime);

    GMMStruct_prime = GMMStruct;

end

GMMStruct_K = GMMStruct;
