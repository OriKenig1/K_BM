
function [GMMStruct,h_opt] = K_BM(x,M,I,GMMStruct_0,MaxIter,tol)

% This function performs GMM parameter estimation using the K-BM algorithm.

% Usage:

% [GMMStruct,h_opt] = K_BM(x,M,I,GMMStruct_0,MaxIter,tol)

% Inputs:

% x - Nxp obserations array, where N denotes the number of observations and p denotes the dimension.

% I - Selection interval for h_opt

% GMMStruct_0 - Initalized GMM parameters structure:

%                            GMMStruct.Dim: Dimension.

%                            GMMStruct.Order: Order (number of Gaussians).

%                            GMMStruct.Alpha: mixing proportions (Orderx1).

%                            GMMStruct.Means: mean vectors (DimensionxOrder).

%                            GMMStruct.Covars: covariance matrices  (DxDxOrder).

% MaxIter - maximal iterations.

% tol - The lgorithm's tolerance between each iterations

% Outputs:

% GMMStruct - Result structure

% h_opt - The selected optimal bandwidth parameter

% By O. Kenig (orikenig@gmail.com) and K. Todros (todros@post.bgu.ac.il)

    arguments
        x double
        M double
        I double = 1:20
        GMMStruct_0 = KMEANS_initGuess(x,M);
        MaxIter double = 50
        tol double = 10^-13
    end

    num_h = length(I);

    parfor h_idx = 1:num_h
             
            curr_h = I(h_idx);

            GMMStruct = K_BM_GMM(x,curr_h,GMMStruct_0,MaxIter,tol);

            AMISE{h_idx} = calc_V(x,GMMStruct);

            results{h_idx} = GMMStruct;

    end

    [~,h_idx] = min( cell2mat(AMISE) ); 

    GMMStruct = cell2mat(results(h_idx));

    h_opt = I(h_idx);
    
end