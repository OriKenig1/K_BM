function [covars,exit] = calc_covars(x,h,GMMStruct,omega,gamma,rho,lambda)

    exit = true;

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    covars = zeros(p,p,M);

    covars_prev = GMMStruct.Covars;

    covars_h = calc_covars_h(covars_prev,h);

    covars_h_inv = calc_covars_inv(covars_h);

    for m = 1:M
    
        covars_w = calc_covarsw(x,GMMStruct.Means(:,m),omega,gamma(:,m),lambda(m));

        covars_phi = calc_covarsphi(x,GMMStruct.Means(:,m),covars_h(:,:,m));

        covars(:,:,m) = covars_w + ...
            rho(m)./lambda(m) .* covars_prev(:,:,m) * covars_h_inv(:,:,m) * ...
            ( covars_h(:,:,m) - covars_phi ) * ...
            covars_h_inv(:,:,m) * covars_prev(:,:,m); 

        curr_det = det(covars(:,:,m));
        curr_cond = cond(covars(:,:,m));
        pos_eig = all(eig(covars(:,:,m)) > 0);

        if curr_cond > 10^3 || abs(curr_det) > 10^3 || ~ pos_eig

            covars(:,:,m) = covars_prev(:,:,m);

        else

            exit = false;

        end

    end
    
    covars = fix_covars(covars,M);

end

