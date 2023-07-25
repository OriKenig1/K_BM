function y = calc_y(GMMStruct,N)

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    y = zeros(GMMStruct.len,1);
  
    if M > 1
        y(1:M-1) = sum(GMMStruct.phi_grad_w_h,2) ./ N;
    end
    
    for m = 1:M

        y(M+(m-1)*p:M+m*p-1) = GMMStruct.Alpha(m)*sum(GMMStruct.phi_grad_u_h(:,:,m,:),4) ./ N;

        y(M+M*p+(m-1)*p*(p+1)/2:M+M*p+m*p*(p+1)/2-1) = GMMStruct.Alpha(m)*sum(GMMStruct.phi_grad_s_h(:,:,m,:),4)./ N;    

    end

    y = y ./ GMMStruct.u;

    % y := grad_u / u
