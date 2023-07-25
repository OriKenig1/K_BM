function F = calc_F(GMMStruct,N)

    F = zeros(GMMStruct.len);

    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    % Pluggin results in F

    F(1:M-1,1:M-1) = zeros(M-1);

    F(1:M-1,M:(M+M*p)-1) = calc_wuh(N,M,p,GMMStruct.phi_grad_u_h);

    F(M:(M+M*p)-1,1:M-1) = F(1:M-1,M:(M+M*p)-1)';

    F(1:M-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_wsh(N,M,p,GMMStruct.phi_grad_s_h);

    F((M+M*p):(M+M*p+M*p*(p+1)/2)-1,1:M-1) = F(1:M-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1)';

    F(M:(M+M*p)-1,M:(M+M*p)-1) = calc_uuh(N,M,p,GMMStruct.Alpha,GMMStruct.phi_grad_uu_h);

    F(M:(M+M*p)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_ush(N,M,p,GMMStruct.Alpha,GMMStruct.phi_grad_us_h);

    F((M+M*p):(M+M*p+M*p*(p+1)/2)-1,M:(M+M*p)-1) = F(M:(M+M*p)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1)';

    F((M+M*p):(M+M*p+M*p*(p+1)/2)-1,(M+M*p):(M+M*p+M*p*(p+1)/2)-1) = calc_ssh(N,M,p,GMMStruct.Alpha,GMMStruct.phi_grad_ss_h);
 
    % F := grad^2_u