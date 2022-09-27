function means = calc_means(x,h,GMMStruct,rho,lambda,omega,gamma)

    M = GMMStruct.Order;

    means = zeros(GMMStruct.Dim, M);

    covars_inv = calc_covars_inv(GMMStruct.Covars);

    covars_h = calc_covars_h(GMMStruct.Covars,h);

    covars_h_inv = calc_covars_inv(covars_h);

    for m = 1:M

        Gamma = covars_inv(:,:,m) - rho(m)./lambda(m) * covars_h_inv(:,:,m);

        mu_omega = ( sum( omega .* gamma(:,m) .* x) ./ ( lambda(m) + 1e-200 ) )';

        phi = mvnpdf(x,GMMStruct.Means(:,m)', covars_h(:,:,m));

        mu_phi = ( sum( phi .* x) ./ ( sum(phi) + 1e-200 ) )';

        means(:,m) = inv(Gamma) * ( covars_inv(:,:,m) * mu_omega - rho(m)/( lambda(m) + 1e-200 ) * covars_h_inv(:,:,m) * mu_phi );

    end

end

