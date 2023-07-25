function  integral_uu = calc_integral_uu(GMMStruct)

    M = GMMStruct.Order;
    
    p = GMMStruct.Dim;

    integral_uu = zeros(M*p,M*p);
    
   for m = 1:M

        for k = 1:M

            tmp = ( GMMStruct.Covars_mk(:,:,m,k) + ( GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m) )*( GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k) )' );

            integral_uu(p*(m-1)+1:p*m,p*(k-1)+1:p*k) = GMMStruct.Alpha(k) * GMMStruct.Alpha(m) * GMMStruct.Gaussians_mk(m,k)...
                * GMMStruct.Covars_inv(:,:,m) * tmp * GMMStruct.Covars_inv(:,:,k) ;

        end

    end
    