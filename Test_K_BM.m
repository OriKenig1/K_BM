clear 
close all

Dimension = 2;
Order = 3;
N = 300;

MaxIter = 30;
tol = 10^-13;

dim_vec = 1:Dimension;

% -------------------------------------------
% Original PDF

mu_1 = -10*ones(Dimension,1);
mu_2 = 3*ones(Dimension,1);
mu_3 = 13*ones(Dimension,1);
mu_Org = [mu_1,mu_2,mu_3];

Sigma_1 = 0.95 .^ (abs( dim_vec'-dim_vec)) ;
Sigma_2 = eye(Dimension);
Sigma_3 = 0.6 .^ (abs( dim_vec'-dim_vec)) ;
Sigma_Org(:,:,1) = Sigma_1;
Sigma_Org(:,:,2) = Sigma_2;
Sigma_Org(:,:,3) = Sigma_3;

alpha_Org = [0.3,0.3,0.4];

GMMStruct_Org = CreateGMMStruct(Dimension, Order);
GMMStruct_Org.Alpha = alpha_Org;
GMMStruct_Org.Means = mu_Org;
GMMStruct_Org.Covars = Sigma_Org;

% -------------------------------------------

% -------------------------------------------
% Contamination PDF

Order_cont = 2;

mu_1 = -20*ones(Dimension,1);
mu_2 = 20*ones(Dimension,1);
mu_cont = [mu_1,mu_2];

Sigma_1 = 40 * inv( 0.6 .^ (abs( dim_vec'-dim_vec)) );
Sigma_2 = 40 * inv( 0.6 .^ (abs( dim_vec'-dim_vec)) );
Sigma_cont(:,:,1) = Sigma_1;
Sigma_cont(:,:,2) = Sigma_2;

alpha_cont = [0.5,0.5];

GMMStruct_cont = CreateGMMStruct(Dimension, Order_cont);
GMMStruct_cont.Alpha = alpha_cont;
GMMStruct_cont.Means = mu_cont;
GMMStruct_cont.Covars = Sigma_cont;

epsilon = 0.1;

% -------------------------------------------
% Initialize 

x = Generate_Obs(GMMStruct_Org,GMMStruct_cont,N,epsilon);

GMMStruct_0 = KMEANS_initGuess(x,Order);

I = 1:1:30;

% ----------------------------------------- %
% Run K-BM %

[GMMStruct,h] = K_BM(x,Order,I,GMMStruct_0);
MISE = MISE(GMMStruct_Org,GMMStruct);
display(MISE,'K-BM MISE');
display(h,'h');
