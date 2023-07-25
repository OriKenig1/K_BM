function G = calc_G(n,GMMStruct)

    G = zeros(GMMStruct.len);

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    % Pluggin results in G

    G(1:M-1,1:M-1) = zeros(M-1);

    G(1:M-1,M:(M+M*p)-1) = calc_wu(M,p,GMMStruct.phi_grad_u(:,:,:,n));

    G(M:(M+M*p)-1,1:M-1) = G(1:M-1,M:(M+M*p)-1)';

    G(1:M-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_ws(M,p,GMMStruct.phi_grad_s(:,:,:,n));

    G((M+M*p):(M+M*p+M*p*(p+1)/2)-1,1:M-1) = G(1:M-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1)';

    G(M:(M+M*p)-1,M:(M+M*p)-1) = calc_uu(M,p,GMMStruct.Alpha,GMMStruct.phi_grad_uu(:,:,:,n));

    G(M:(M+M*p)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_us(M,p,GMMStruct.Alpha,GMMStruct.phi_grad_us(:,:,:,n));

    G((M+M*p):(M+M*p+M*p*(p+1)/2)-1,M:(M+M*p)-1) = G(M:(M+M*p)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1)';

    G((M+M*p):(M+M*p+M*p*(p+1)/2)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_ss(M,p,GMMStruct.Alpha,GMMStruct.phi_grad_ss(:,:,:,n));
