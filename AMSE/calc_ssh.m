function ssh = calc_ssh(N,M,p,alpha,phi_grad_ss_h)

    ssh = zeros(M*p*(p+1)/2,M*p*(p+1)/2);

    for m = 1:M

        for k = 1:M

            F = (m==k)*alpha(k)*sum(phi_grad_ss_h(:,:,k,:),4)' ./ N;

            ssh(p*(p+1)/2*(m-1)+1:p*(p+1)/2*m,p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = F;

        end
        
    end