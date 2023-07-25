function KBIC = calc_KBIC(x,GMMStruct,h,M)
    
    [N,p] = size(x);

    w = calc_w(x,h);

    Gaussians = zeros(M,1);
    Gaussians_h = zeros(M,1);

    f = zeros(N,1);
    f_h = zeros(N,1);

    for l = 1:M

        GMMStruct.Covars_h(:,:,l) = GMMStruct.Covars(:,:,l) + h^2*eye(p);

    end

    for n = 1:N

        for m = 1:M

            Gaussians(m) = mvnpdf(x(n,:),GMMStruct.Means(:,m)',GMMStruct.Covars(:,:,m));
            Gaussians_h(m) = mvnpdf(x(n,:),GMMStruct.Means(:,m)',GMMStruct.Covars_h(:,:,m));

        end

        f(n) = GMMStruct.Alpha'*Gaussians;
        f_h(n) = GMMStruct.Alpha'*Gaussians_h;

    end

    m = GMMStruct.len;

    u = sum(f_h) ./ N;

    J = sum(w.*log(f)) - log(u);

    KBIC = -J + 0.5*m/N*log(N);
