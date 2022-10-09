function fixed_covars = fix_covars(covars,M)

    fixed_covars = covars;

    for m = 1:M

        fixed_covars(:,:,m) = (fixed_covars(:,:,m) + fixed_covars(:,:,m).')./2;

    end

