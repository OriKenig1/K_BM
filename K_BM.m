
function [GMMStruct,h] = K_BM(x,M,I,GMMStruct_0,MaxIter,tol)

    % K_BM, runs the K-BM algorithm with data x. 
    % Selection of the bandiwdth parameter is done via the L2 method.
    % M is the number of clusters, I is the search interval of h
    % GMMStrcuct_ 0 is the initialized value of the vector parameters,
    % default is calculated via the k-means method.
    % The function returns the estimated GMM and the selectd bandwidth
    % parameter

    arguments
        x double
        M double
        I double = 1:20
        GMMStruct_0 = KMEANS_initGuess(x,M);
        MaxIter double = 30
        tol double = 10^-13
    end

    num_h = length(I);

    parfor h_idx = 1:num_h
             
            curr_h = I(h_idx);

            GMMStruct = K_BM_GMM(x,M,curr_h,GMMStruct_0,MaxIter,tol);

            AMISE{h_idx} = calc_V(x,GMMStruct);

            results{h_idx} = GMMStruct;

    end

    [~,h_idx] = min( cell2mat(AMISE) ); 

    GMMStruct = cell2mat(results(h_idx));

    h = I(h_idx);
    
end