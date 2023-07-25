function  integral_wu = calc_integral_wu(GMMStruct)

    M = GMMStruct.Order;
    
    p = GMMStruct.Dim;

    integral_wu = zeros(M-1,M*p);
    
   for m = 2:M

        for k = 1:M

            el1 = GMMStruct.Gaussians_mk(m,k) * ( GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k) )';

            el2 = GMMStruct.Gaussians_mk(1,k) * ( GMMStruct.Means_mk(:,:,1,k) - GMMStruct.Means(:,k) )';

            integral_wu(m-1,p*(k-1)+1:p*k) = GMMStruct.Alpha(k) * ( el1 - el2 ) * GMMStruct.Covars_inv(:,:,k);

        end

    end
    