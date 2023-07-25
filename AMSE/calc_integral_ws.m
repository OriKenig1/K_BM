function  integral_ws = calc_integral_ws(GMMStruct)

    M = GMMStruct.Order;
    
    p = GMMStruct.Dim;

    integral_ws = zeros(M-1,M*p*(p+1)/2);
   
    L = get_L(p);

    for m = 2:M

        for k = 1:M

            el1 = GMMStruct.Covars_inv(:,:,k) * ( GMMStruct.Gaussians_mk(m,k) - GMMStruct.Gaussians_mk(1,k) );

            el2 = GMMStruct.Gaussians_mk(m,k) * (GMMStruct.Covars_mk(:,:,m,k) + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))');

            el3 = GMMStruct.Gaussians_mk(1,k) * (GMMStruct.Covars_mk(:,:,1,k) + (GMMStruct.Means_mk(:,:,1,k) - GMMStruct.Means(:,k))*(GMMStruct.Means_mk(:,:,1,k) - GMMStruct.Means(:,k))');

            vec = el1 - GMMStruct.Covars_inv(:,:,k) * ( el2 - el3 ) * GMMStruct.Covars_inv(:,:,k);

            integral_ws(m-1,p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = -0.5 * GMMStruct.Alpha(k) * vec(:)' * L';

        end

    end
    