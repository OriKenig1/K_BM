function covars_h = calc_covars_h(covars,h)

    [~,p,M] = size(covars);

    covars_h = zeros(p,p,M);

    for m = 1:M
        covars_h(:,:,m) = covars(:,:,m) + h^2 * eye(p);
    end

end