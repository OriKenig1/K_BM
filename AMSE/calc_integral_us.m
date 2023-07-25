function  integral_us = calc_integral_us(GMMStruct)

    M = GMMStruct.Order;
    
    p = GMMStruct.Dim;

    integral_us = zeros(M*p,M*p*(p+1)/2);
    
   for m = 1:M

        for k = 1:M

            cov_inv_k = GMMStruct.Covars_inv(:,:,k);
         
            L = get_L(p);

            el1 = ( GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m) )  * cov_inv_k(:)';
            
            el2 = GMMStruct.Covars_mk(:,:,m,k)*kron(GMMStruct.Covars_inv(:,:,k), (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))'*GMMStruct.Covars_inv(:,:,k));

            el3 = GMMStruct.Covars_mk(:,:,m,k)*kron((GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))'*GMMStruct.Covars_inv(:,:,k), GMMStruct.Covars_inv(:,:,k));

            tmp1 = GMMStruct.Covars_inv(:,:,k)*GMMStruct.Covars_mk(:,:,m,k)*GMMStruct.Covars_inv(:,:,k);

            el4 = (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))*tmp1(:)';

            tmp2 = GMMStruct.Covars_inv(:,:,k)*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))'*GMMStruct.Covars_inv(:,:,k);

            el5 = (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))*tmp2(:)';

            integral_us(p*(m-1)+1:p*m,p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = -0.5*GMMStruct.Alpha(m) * GMMStruct.Alpha(k) * GMMStruct.Gaussians_mk(m,k) * GMMStruct.Covars_inv(:,:,m)*...
                ( el1 - ( el2 + el3 + el4 + el5 ) )*L';

        end

    end
    