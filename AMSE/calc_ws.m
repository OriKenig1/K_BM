function ws = calc_ws(M,p,phi_grad_s)

    ws = zeros(M-1,M*p*(p+1)/2);

    for m = 2:M

        for k = 1:M

            G = (m==k)*phi_grad_s(:,:,k)';

            if k == 1
                G = - phi_grad_s(:,:,1)';
            end

            ws((m-1),p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = G;

        end

    end