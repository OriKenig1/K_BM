function  GMMStruct = calculate_derivatives(all_x,GMMStruct)

    N = length(all_x);
    
    M = GMMStruct.Order;

    p = GMMStruct.Dim;

    L = get_L(p);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    GMMStruct.phi_grad_u = zeros(p,1,M,N);
    GMMStruct.phi_grad_u_h = zeros(p,1,M,N);
    GMMStruct.phi_grad_s = zeros(p*(p+1)/2,1,M,N);
    GMMStruct.phi_grad_s_h = zeros(p*(p+1)/2,1,M,N);
    GMMStruct.phi_grad_w = zeros(M-1,N);
    GMMStruct.phi_grad_w_h = zeros(M-1,N);
    GMMStruct.phi_grad_uu = zeros(p,p,M,N);
    GMMStruct.phi_grad_uu_h = zeros(p,p,M,N);
    GMMStruct.phi_grad_us = zeros(p,p*(p+1)/2,M,N);
    GMMStruct.phi_grad_us_h = zeros(p,p*(p+1)/2,M,N);
    GMMStruct.phi_grad_ss = zeros(p*(p+1)/2,p*(p+1)/2,M,N);
    GMMStruct.phi_grad_ss_h = zeros(p*(p+1)/2,p*(p+1)/2,M,N);
    GMMStruct.q = zeros(GMMStruct.len,N);
    GMMStruct.b = zeros(GMMStruct.len,N);

    for n = 1:N

        x = all_x(n,:)';
        
        % Calculating derivatives of gaussians w.r.t. mean and covariance
        for m = 1:M
            GMMStruct.phi_grad_u(:,:,m,n) = GMMStruct.Gaussians(n,m)*GMMStruct.Covars_inv(:,:,m)*(x-GMMStruct.Means(:,m));
            GMMStruct.phi_grad_u_h(:,:,m,n) = GMMStruct.Gaussians_h(n,m)*GMMStruct.Covars_h_inv(:,:,m)*(x-GMMStruct.Means(:,m));
            tmp = GMMStruct.Covars_inv(:,:,m) - calc_abba(GMMStruct.Covars_inv(:,:,m),x-GMMStruct.Means(:,m));
            GMMStruct.phi_grad_s(:,:,m,n) = -0.5*GMMStruct.Gaussians(n,m)*L*tmp(:);
            tmp = GMMStruct.Covars_h_inv(:,:,m) - calc_abba(GMMStruct.Covars_h_inv(:,:,m),x-GMMStruct.Means(:,m));
            GMMStruct.phi_grad_s_h(:,:,m,n) = -0.5*GMMStruct.Gaussians_h(n,m)*L*tmp(:);
        end

        % Calculating derivatives of f w.r.t. omega
        for m = 2:M
            GMMStruct.phi_grad_w(m-1,n) = GMMStruct.Gaussians(n,m) - GMMStruct.Gaussians(n,1);
            GMMStruct.phi_grad_w_h(m-1,n) = GMMStruct.Gaussians_h(n,m) - GMMStruct.Gaussians_h(n,1);
        end

        % Calculating second derivatives of gaussians w.r.t. mean and sigma
        for m = 1:M
            %%%%%%%%%%%%%%%%%%
            GMMStruct.phi_grad_uu(:,:,m,n) = GMMStruct.Gaussians(n,m)*...
                (calc_abba(GMMStruct.Covars_inv(:,:,m),x-GMMStruct.Means(:,m)) - GMMStruct.Covars_inv(:,:,m));
            GMMStruct.phi_grad_uu_h(:,:,m,n) = GMMStruct.Gaussians_h(n,m)*...
                (calc_abba(GMMStruct.Covars_h_inv(:,:,m),x-GMMStruct.Means(:,m)) - GMMStruct.Covars_h_inv(:,:,m));
            %%%%%%%%%%%%%%%%%%
            tmp = GMMStruct.Covars_inv(:,:,m) - calc_abba(GMMStruct.Covars_inv(:,:,m),x-GMMStruct.Means(:,m));
            GMMStruct.phi_grad_us(:,:,m,n) = -0.5*GMMStruct.Gaussians(n,m)*GMMStruct.Covars_inv(:,:,m)...
                *( (x-GMMStruct.Means(:,m))*tmp(:)' + 2*kron( (x-GMMStruct.Means(:,m))'*GMMStruct.Covars_inv(:,:,m), eye(p) ) ) * L';
            tmp = GMMStruct.Covars_h_inv(:,:,m) - calc_abba(GMMStruct.Covars_h_inv(:,:,m),x-GMMStruct.Means(:,m));
            GMMStruct.phi_grad_us_h(:,:,m,n) = -0.5*GMMStruct.Gaussians_h(n,m)*GMMStruct.Covars_h_inv(:,:,m)...
                *( (x-GMMStruct.Means(:,m))*tmp(:)' + 2*kron( (x-GMMStruct.Means(:,m))'*GMMStruct.Covars_h_inv(:,:,m), eye(p) ) ) * L';
            %%%%%%%%%%%%%%%%%%
            tmp = GMMStruct.Covars_inv(:,:,m) - calc_abba(GMMStruct.Covars_inv(:,:,m),x-GMMStruct.Means(:,m));
            tmp2 = calc_abba(GMMStruct.Covars_inv(:,:,m),x-GMMStruct.Means(:,m));
            GMMStruct.phi_grad_ss(:,:,m,n) = 0.25*GMMStruct.Gaussians(n,m)*L*...
                    ( tmp(:)*tmp(:)' + 2*(kron(GMMStruct.Covars_inv(:,:,m),GMMStruct.Covars_inv(:,:,m)) - kron(GMMStruct.Covars_inv(:,:,m),tmp2)...
                        - kron(tmp2,GMMStruct.Covars_inv(:,:,m))) )*L';
            tmp = GMMStruct.Covars_h_inv(:,:,m) - calc_abba(GMMStruct.Covars_h_inv(:,:,m),x-GMMStruct.Means(:,m));
            tmp2 = calc_abba(GMMStruct.Covars_h_inv(:,:,m),x-GMMStruct.Means(:,m));
            GMMStruct.phi_grad_ss_h(:,:,m,n) = 0.25*GMMStruct.Gaussians_h(n,m)*L*...
                    ( tmp(:)*tmp(:)' + 2*(kron(GMMStruct.Covars_h_inv(:,:,m),GMMStruct.Covars_h_inv(:,:,m)) - kron(GMMStruct.Covars_h_inv(:,:,m),tmp2)...
                        - kron(tmp2,GMMStruct.Covars_h_inv(:,:,m))) )*L';
        end

        GMMStruct.q(:,n) = calc_q(n,GMMStruct);
        GMMStruct.b(:,n) = calc_b(n,GMMStruct);

    end

    GMMStruct.y = calc_y(GMMStruct,N);

end

function abba = calc_abba(a,b)

    tmp = a*b;
    abba = tmp*tmp';

end
