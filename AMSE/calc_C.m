function C = calc_C(x,h,GMMStruct)

    N = length(x);

    w = calc_w(x,h);

    % Derivatives of f

    H = zeros(GMMStruct.len,GMMStruct.len,N);
    
    for n = 1:N

        H(:,:,n) = calc_G(n,GMMStruct) ./ (GMMStruct.f(n)+1e-10) - GMMStruct.q(:,n)*GMMStruct.q(:,n)';
        
    end

    % Derivatives of u

    hess_logu = calc_F(GMMStruct,N) ./ GMMStruct.u - GMMStruct.y*GMMStruct.y';

    % Plug into C

    C = sum(H.*reshape(w,[1,1,N]),3) - hess_logu;