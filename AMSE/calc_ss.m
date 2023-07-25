function ss = calc_ss(M,p,alpha,phi_grad_ss)

    ss = zeros(M*p*(p+1)/2,M*p*(p+1)/2);

    for m = 1:M

        for k = 1:M

            G = (m==k)*alpha(k)*phi_grad_ss(:,:,k);

            ss(p*(p+1)/2*(m-1)+1:p*(p+1)/2*m,p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = G;
            
        end

    end