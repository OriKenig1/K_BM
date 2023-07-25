function uu = calc_uu(M,p,alpha,phi_grad_uu)

    uu = zeros(M*p,M*p);

    for m = 1:M

        for k = 1:M

            G = (m==k)*alpha(k)*phi_grad_uu(:,:,k);

            uu(p*(m-1)+1:p*m,p*(k-1)+1:p*k) = G;

        end

    end