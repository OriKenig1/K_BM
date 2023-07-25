function  integral_ww = calc_integral_ww(GMMStruct)

    M = GMMStruct.Order;

    integral_ww = zeros(M-1);
    
   for m = 2:M

        for k = 2:M

            integral_ww(m-1,k-1) = GMMStruct.Gaussians_mk(m,k) - GMMStruct.Gaussians_mk(m,1)...
                - GMMStruct.Gaussians_mk(k,1) + GMMStruct.Gaussians_mk(1,1);

        end

    end
    