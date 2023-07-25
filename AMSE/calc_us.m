function us = calc_us(M,p,alpha,phi_grad_us)

    us = zeros(M*p,M*p*(p+1)/2);

    for m = 1:M

        for k = 1:M

            G = (m==k)*alpha(k)*phi_grad_us(:,:,k);

            us(p*(m-1)+1:p*m,p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = G;

        end

    end