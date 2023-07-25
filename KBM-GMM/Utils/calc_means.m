function means = calc_means(x,h,GMMStruct,rho,lambda,weights,gamma)

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    means = zeros(GMMStruct.Dim, M);

    covars_inv = calc_covars_inv(GMMStruct.Covars);

    covars_h = calc_covars_h(GMMStruct.Covars,h);

    covars_h_inv = calc_covars_inv(covars_h);

    for m = 1:M

        Gamma = covars_inv(:,:,m) - rho(m)/( lambda(m) + 1e-20 ) * covars_h_inv(:,:,m);

        invGamma = inv(Gamma + eye(p)*1e-20);

        mu_omega = ( sum( weights .* gamma(:,m) .* x) ./ ( lambda(m) + 1e-20 ) )';

        phi = mvnpdf(x,GMMStruct.Means(:,m)', covars_h(:,:,m));

        mu_phi = ( sum( phi .* x) ./ ( sum(phi) + 1e-20 ) )';

        means(:,m) = invGamma * ( covars_inv(:,:,m) * mu_omega - rho(m)/( lambda(m) + 1e-200 ) * covars_h_inv(:,:,m) * mu_phi );

    end

end

