function [K_GMMStruct,idx_l_opt,idx_h_opt] = COMPLETE_KBM(x,h,L,MaxIter,tol)

% Running the K-BM on all avilable bandwidths and model orders
fprintf("Running the K-BM on all model orders:");

for l_idx = 1:length(L) 

    curr_L = L(l_idx);
    
    GMMStruct_0 = KMEANS_initGuess(x,curr_L);

    parfor h_idx = 1:length(h)

        curr_h = h(h_idx);
        
        GMMStruct(h_idx,l_idx) = K_BM_GMM(x,curr_h,GMMStruct_0,MaxIter,tol);

    end

    fprintf(" %d...", curr_L);

end

fprintf("\nCalculating the best model order for each bandwidth parameter...\n");

optimal_Order = zeros(length(h),1);

for h_idx = 1:length(h)

    curr_h = h(h_idx);

    BIC = zeros(length(L),1);

    parfor l_idx = 1:length(L)

        curr_L = L(l_idx);
        
        curr_GMMStruct = GMMStruct(h_idx,l_idx);

        BIC(l_idx) = calc_KBIC(x,curr_GMMStruct,curr_h,curr_L);

    end

    [~, bic_opt_idx] = min(BIC);

    optimal_Order(h_idx) = bic_opt_idx;

end

fprintf("Finding the best (bandwidth,order) w.r.t. the MISE\n");
Vhat = zeros(length(h),1);

parfor h_idx = 1:length(h)

    curr_h = h(h_idx);

    curr_l_idx = optimal_Order(h_idx);

    Vhat(h_idx)= calc_AMSE(x, curr_h, GMMStruct(h_idx,curr_l_idx));

end

[~, idx_h_opt] = min(Vhat);

% Final Results %

K_GMMStruct = GMMStruct(idx_h_opt,optimal_Order(idx_h_opt));

idx_l_opt = optimal_Order(idx_h_opt);







