function wsh = calc_wsh(N,M,p,phi_grad_s_h)

    wsh = zeros(M-1,M*p*(p+1)/2);

    for m = 2:M

        for k = 1:M

            F = (m==k)*sum(phi_grad_s_h(:,:,k,:),4)' ./ N;

            if k == 1
                F = - sum(phi_grad_s_h(:,:,1,:),4)' ./ N;
            end

            wsh((m-1),p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = F;


        end

    end