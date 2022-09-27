function fixed_covars = fix_covars(covars,M)

    %[p,~,~] = size(covars);

    fixed_covars = covars;

    for m = 1:M

        fixed_covars(:,:,m) = (fixed_covars(:,:,m) + fixed_covars(:,:,m).')./2;

        %T = real(sqrtm(fixed_covars(:,:,m)));

        %fixed_covars(:,:,m) = T*T.';

    end

