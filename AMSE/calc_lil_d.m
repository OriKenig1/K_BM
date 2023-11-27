function d = calc_lil_d(x,h,GMMStruct,c,sum_g_tild)

    N = length(x);

    p = GMMStruct.Dim;

    d = zeros(GMMStruct.len,N);

    K_h0 = mvnpdf(zeros(1,p),zeros(1,p),h^2*eye(p));

    phi_0 = K_h0 ./ ( N * sum_g_tild );

    for n = 1:N

        K_h = mvnpdf(x - x(n,:),zeros(1,p),h^2*eye(p));

        phi = K_h ./ ( N * sum_g_tild );
   
        d(:,n) = (N-1) * sum( phi'.*c - c(:,n)*phi_0./N , 2);

    end
