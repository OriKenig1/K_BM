function b = calc_b(n,GMMStruct)

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    b = zeros(GMMStruct.len,1);
  
    if M > 1
        b(1:M-1) = GMMStruct.phi_grad_w_h(:,n);
    end
    
    for m = 1:M

        b(M+(m-1)*p:M+m*p-1) = GMMStruct.Alpha(m).*GMMStruct.phi_grad_u_h(:,:,m,n);

        b(M+M*p+(m-1)*p*(p+1)/2:M+M*p+m*p*(p+1)/2-1) = GMMStruct.Alpha(m).*GMMStruct.phi_grad_s_h(:,:,m,n);

    end

    b = b ./ (GMMStruct.f_h(n));
    
    % b = grad_fh ./ fh