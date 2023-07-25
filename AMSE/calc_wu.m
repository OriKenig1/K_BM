function wu = calc_wu(M,p,phi_grad_u)

    wu = zeros(M-1,M*p);

    for m = 2:M

        for k = 1:M

            G = (m==k)*phi_grad_u(:,:,k)';

            if k == 1
                G = - phi_grad_u(:,:,1)';
            end

            wu((m-1),p*(k-1)+1:p*k) = G;

        end

    end