function  Integral = calc_Integral(GMMStruct)

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    % Calculate common mean and variance
    for m = 1:M
        for k = 1:M

            GMMStruct.Gaussians_mk(m,k) = mvnpdf(GMMStruct.Means(:,m)',GMMStruct.Means(:,k)',GMMStruct.Covars(:,:,m)+GMMStruct.Covars(:,:,k));

            GMMStruct.Covars_mk(:,:,m,k) = inv( GMMStruct.Covars_inv(:,:,m) + GMMStruct.Covars_inv(:,:,k) + 1e-3*eye(p) );

            GMMStruct.Means_mk(:,:,m,k) = GMMStruct.Covars_mk(:,:,m,k) * ...
                ( GMMStruct.Covars_inv(:,:,m) * GMMStruct.Means(:,m) + GMMStruct.Covars_inv(:,:,k) * GMMStruct.Means(:,k) );

        end
    end

    % Fix floating point errors for mean and covariance of same index    
    for k = 1:M
        GMMStruct.Gaussians_mk(k,k) = 1 ./ ( sqrt( (2*pi)^(-p)*2 * det(GMMStruct.Covars(:,:,k)) ) + 1e-10 );
        GMMStruct.Covars_mk(:,:,k,k) = 0.5 * GMMStruct.Covars(:,:,k);
        GMMStruct.Means_mk(:,:,k,k) = GMMStruct.Means(:,k);          
    end

    Integral = zeros(GMMStruct.len);
    
    % Start filling Integral

    Integral(1:M-1,1:M-1) = calc_integral_ww(GMMStruct);

    Integral(1:M-1,M:(M+M*p)-1) = calc_integral_wu(GMMStruct);

    Integral(M:(M+M*p)-1,1:M-1) = Integral(1:M-1,M:(M+M*p)-1)';

    Integral(1:M-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_integral_ws(GMMStruct);

    Integral((M+M*p):(M+M*p+M*p*(p+1)/2)-1,1:M-1) = Integral(1:M-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1)';

    Integral(M:(M+M*p)-1,M:(M+M*p)-1) = calc_integral_uu(GMMStruct);

    Integral(M:(M+M*p)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_integral_us(GMMStruct);

    Integral((M+M*p):(M+M*p+M*p*(p+1)/2)-1,M:(M+M*p)-1) = Integral(M:(M+M*p)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1)';

    Integral((M+M*p):(M+M*p+M*p*(p+1)/2)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_integral_ss(GMMStruct);

    Integral = ( Integral + Integral' ) ./ 2;
