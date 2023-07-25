function wuh = calc_wuh(N,M,p,phi_grad_u_h)

    wuh = zeros(M-1,M*p);

    for m = 2:M

        for k = 1:M

            F = (m==k)*sum(phi_grad_u_h(:,:,k,:),4)' ./ N;

            if k == 1
                F = - sum(phi_grad_u_h(:,:,1,:),4)' ./ N;
            end

            wuh((m-1),p*(k-1)+1:p*k) = F;
        end

    end