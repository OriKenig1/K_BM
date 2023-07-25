function  integral_ss = calc_integral_ss(GMMStruct)

   M = GMMStruct.Order;
    
   p = GMMStruct.Dim;

   integral_ss = zeros(M*p*(p+1)/2);
    
   e=@(a) [zeros(a-1,1);1;zeros(p-a,1)];

   L = get_L(p);

   for m = 1:M

        covars_inv_m = GMMStruct.Covars_inv(:,:,m);

        for k = 1:M

            covars_inv_k = GMMStruct.Covars_inv(:,:,k);

            J = covars_inv_m(:) * covars_inv_k(:)';

            tmp1 = GMMStruct.Covars_inv(:,:,k) * ( GMMStruct.Covars_mk(:,:,m,k) + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k)) * (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))' ) * GMMStruct.Covars_inv(:,:,k);

            P = covars_inv_m(:) * tmp1(:)';

            tmp2 = GMMStruct.Covars_inv(:,:,m) * ( GMMStruct.Covars_mk(:,:,m,k) + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m)) * (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))' ) * GMMStruct.Covars_inv(:,:,m);

            Q = tmp2(:) * covars_inv_k(:)';

            R = 0;

            for i = 1:p
                for j = 1:p

                G = GMMStruct.Covars_inv(:,:,m)*e(i)*e(j)'*GMMStruct.Covars_inv(:,:,k);

                Psi = GMMStruct.Covars_mk(:,:,m,k)*(G+G')*GMMStruct.Covars_mk(:,:,m,k) + trace(G*GMMStruct.Covars_mk(:,:,m,k))*GMMStruct.Covars_mk(:,:,m,k) ...
                    + GMMStruct.Covars_mk(:,:,m,k)*G*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))' ...
                    + GMMStruct.Covars_mk(:,:,m,k)*G'*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))' ...
                    + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))'*G*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))*GMMStruct.Covars_mk(:,:,m,k) ...
                    + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))'*trace(G*GMMStruct.Covars_mk(:,:,m,k)) ...
                    + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))'*G'*GMMStruct.Covars_mk(:,:,m,k) ...
                    + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))'*G*GMMStruct.Covars_mk(:,:,m,k) ...
                    + (GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,m))'*G*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))*(GMMStruct.Means_mk(:,:,m,k) - GMMStruct.Means(:,k))';

                R = R + kron(e(i),eye(p))*GMMStruct.Covars_inv(:,:,m)*Psi*GMMStruct.Covars_inv(:,:,k)*kron(e(j)',eye(p));

                end
            end

            T = J - P - Q + R;

            integral_ss(p*(p+1)/2*(m-1)+1:p*(p+1)/2*m,p*(p+1)/2*(k-1)+1:p*(p+1)/2*k) = 0.25 * GMMStruct.Alpha(k) * GMMStruct.Alpha(m) * GMMStruct.Gaussians_mk(m,k) * L * T * L';
        end

    end
    