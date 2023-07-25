function GMMStruct = CreateGMMStruct(Dimension, Order)



% This function creates a structure with  GMM parameters.



% Usage:

% GMMStruct = CreateGMMStruct(Dimension, Order)



% Inputs:

% Dimenstion - dimenstion.

% Order - model order (numer of Gaussians).



% Outputs:

% GMMStruct - GMM parameters structure:

%                            GMMStruct.Dim: Dimension.

%                            GMMStruct.Order: Order (number of Gaussians).

%                            GMMStruct.Alpha: mixing proportions (Orderx1).

%                            GMMStruct.Means: mean vectors (DimensionxOrder).

%                            GMMStruct.Covars: covariance matrices  (DxDxOrder).



% Version 1.0, Mai 2006.

% By Koby Todros  todros@widemed.com





GMMStruct.Dim = Dimension;

GMMStruct.Order = Order;

GMMStruct.Alpha = ones(Order,1)./Order;

for i = 1:Order

    GMMStruct.Means(:,i) = randn(1)*ones(Dimension,1);
    
    GMMStruct.Covars(:,:,i) = eye(Dimension,Dimension);

end

GMMStruct.len = (Order - 1) + Order*Dimension + Order*(Dimension*(Dimension+1)/2);


