function q = calc_q(n,GMMStruct)

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    q = zeros(GMMStruct.len,1);
  
    if M > 1
        q(1:M-1) = GMMStruct.phi_grad_w(:,n);
    end
    
    for m = 1:M

        q(M+(m-1)*p:M+m*p-1) = GMMStruct.Alpha(m).*GMMStruct.phi_grad_u(:,:,m,n);

        q(M+M*p+(m-1)*p*(p+1)/2:M+M*p+m*p*(p+1)/2-1) = GMMStruct.Alpha(m).*GMMStruct.phi_grad_s(:,:,m,n);

    end

    q = q ./ (GMMStruct.f(n));
