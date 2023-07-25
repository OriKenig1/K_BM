function ush = calc_ush(N,M,p,alpha,phi_grad_us_h)

    ush = zeros(M*p,M*p*(p+1)/2);

    for m = 1:M

        for k = 1:M

            F = (m==k)*alpha(k)*sum(phi_grad_us_h(:,:,k,:),4) ./ N;

            ush(p*(m-1)+1:p*m,p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = F;

        end

    end