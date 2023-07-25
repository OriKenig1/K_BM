function AMSE = calc_AMSE(x,h,GMMStruct)

    N = length(x);

    M = GMMStruct.Order;

    p = GMMStruct.Dim;
    
    % Prepare variables 

    for m = 1:M

        GMMStruct.Covars_h(:,:,m) = GMMStruct.Covars(:,:,m) + h^2*eye(p);
        GMMStruct.Covars_inv(:,:,m) = inv( GMMStruct.Covars(:,:,m) + 1e-3*eye(p) );
        GMMStruct.Covars_h_inv(:,:,m) = inv( GMMStruct.Covars_h(:,:,m) + 1e-3*eye(p) );

    end


    for n = 1:N

        for m = 1:M

            GMMStruct.Gaussians(n,m) = mvnpdf(x(n,:),GMMStruct.Means(:,m)',GMMStruct.Covars(:,:,m));
            GMMStruct.Gaussians_h(n,m) = mvnpdf(x(n,:),GMMStruct.Means(:,m)',GMMStruct.Covars_h(:,:,m));

        end

    end

    % Fix for too low gaussian pdf values
    min_not_zero = min(GMMStruct.Gaussians(GMMStruct.Gaussians>0));
    GMMStruct.Gaussians(GMMStruct.Gaussians==0) = 1e-40;

    min_not_zero_h = min(GMMStruct.Gaussians_h(GMMStruct.Gaussians_h>0));
    GMMStruct.Gaussians_h(GMMStruct.Gaussians_h==0) = 1e-40;

    for n = 1:N

        GMMStruct.f(n) = sum(GMMStruct.Alpha'.*GMMStruct.Gaussians(n,:));
        GMMStruct.f_h(n) = sum(GMMStruct.Alpha'.*GMMStruct.Gaussians_h(n,:));

    end


    GMMStruct.u = sum(GMMStruct.f_h) ./ N;

    % Start calculations

    [R,~,~] = calc_R(x,h,GMMStruct);

    Integral = calc_Integral(GMMStruct);

    AMSE = trace(R * Integral) + calc_V(x,GMMStruct);
    