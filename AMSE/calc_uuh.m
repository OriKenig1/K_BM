function uuh = calc_uuh(N,M,p,alpha,phi_grad_uu_h)

    uuh = zeros(M*p,M*p);

    for m = 1:M

        for k = 1:M

            F = (m==k)*alpha(k)*sum(phi_grad_uu_h(:,:,k,:),4)' ./ N;

            uuh(p*(m-1)+1:p*m,p*(k-1)+1:p*k) = F;

        end

    end