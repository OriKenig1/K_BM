function covars_h_inv = calc_covars_inv(covars)

    [~,p,M] = size(covars);

    covars_h_inv = zeros(p,p,M);

    for m = 1:M
        covars_h_inv(:,:,m) = inv( covars(:,:,m) + 10^-20 * eye(p) );
    end

end