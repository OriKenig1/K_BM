function covars = calc_covars(x,h,GMMStruct,weights,gamma,rho,lambda)

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    covars = zeros(p,p,M);

    covars_prev = GMMStruct.Covars;

    covars_h = calc_covars_h(covars_prev,h);

    covars_h_inv = calc_covars_inv(covars_h);

    for m = 1:M
    
        covars_w = calc_covarsw(x,GMMStruct.Means(:,m),weights,gamma(:,m),lambda(m));

        covars_phi = calc_covarsphi(x,GMMStruct.Means(:,m),covars_h(:,:,m));
        
        % Update the eigenvectors matrix

        covars(:,:,m) = covars_w + ...
            rho(m)./lambda(m) .* covars_prev(:,:,m) * covars_h_inv(:,:,m) * ...
            ( covars_h(:,:,m) - covars_phi ) * ...
            covars_h_inv(:,:,m) * covars_prev(:,:,m) + 1e-6*eye(p); 

        covars(:,:,m) = ( covars(:,:,m) + covars(:,:,m)' ) ./ 2;
        
        [VEC,EIG] = eig(covars(:,:,m));

        % Update the eigenvalues matrix

        lambda_max = max( std(x).^2 );

        for q = 1:p
            EIG(q,q) = min(max(EIG(q,q),lambda_max/10^5),lambda_max);
        end

        % Reconstruct

        covars(:,:,m) = VEC*EIG*VEC';

        covars(:,:,m) = ( covars(:,:,m) + covars(:,:,m)' ) ./ 2;

    end
    
end

