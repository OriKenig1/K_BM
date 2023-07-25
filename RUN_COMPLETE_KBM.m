clear 
clc
close all

Dimension = 6;
N = 300;

MaxIter = 50;
tol = 10^-6;

epsilon = 0.1;

h = linspace(2,8,20); % Bandwidth parameters
L = 3; % Model orders

% -------------------------------------------
% Original PDF

Order = 3;

dim_vec = 1:Dimension;

v1 = ( (-1).^(dim_vec-1) )';

mu_1 = 2*ones(1,Dimension)';
mu_2 = 10*ones(1,Dimension)';
mu_3 = mu_1 + v1*6;
mu_Org = [mu_1,mu_2,mu_3];

Sigma_1 = 5 * 0.4 .^ (abs( dim_vec'-dim_vec)) ;
Sigma_2 = 3 * eye(Dimension);
Sigma_3 = 3 * 0.6 .^ (abs( dim_vec'-dim_vec)) ;
Sigma_Org(:,:,1) = Sigma_1;
Sigma_Org(:,:,2) = Sigma_2;
Sigma_Org(:,:,3) = Sigma_3;

alpha_Org = [0.3,0.3,0.4]';

GMMStruct_Org = CreateGMMStruct(Dimension, Order);
GMMStruct_Org.Alpha = alpha_Org;
GMMStruct_Org.Means = mu_Org;
GMMStruct_Org.Covars = Sigma_Org;

% -------------------------------------------

% -------------------------------------------
% Contamination PDF

Order_cont = 3;

mu_1 = 10*ones(1,Dimension)' - v1*10;
mu_2 = -20*ones(1,Dimension)' + v1*5;
mu_3 = 5*ones(1,Dimension)' + v1*3.5;
mu_cont = [mu_1,mu_2,mu_3];

Sigma_1 = 10 * 0.9 .^ (abs( dim_vec'-dim_vec));
Sigma_2 = 10 * 0.9 .^ (abs( dim_vec'-dim_vec));
Sigma_3 = 1*eye(Dimension);
Sigma_cont(:,:,1) = Sigma_1;
Sigma_cont(:,:,2) = Sigma_2;
Sigma_cont(:,:,3) = Sigma_3;

alpha_cont = [0.4,0.4,0.2]';

GMMStruct_cont = CreateGMMStruct(Dimension, Order_cont);
GMMStruct_cont.Alpha = alpha_cont;
GMMStruct_cont.Means = mu_cont;
GMMStruct_cont.Covars = Sigma_cont;

% -------------------------------------------

% run K-BM %

[x,x_org,x_cont] = Generate_Obs(GMMStruct_Org,GMMStruct_cont,N,epsilon);

[K_GMMStruct,l_opt_idx,h_opt_idx] = COMPLETE_KBM(x,h,L,MaxIter,tol);

% plot results %

best_Order = L(l_opt_idx);
best_Bandwidth = h(h_opt_idx);

plot_3D(K_GMMStruct,x_org,x_cont,sprintf("%d Clusters ( h = %.1d )",best_Order,best_Bandwidth),'black');




if K_GMMStruct.Order == GMMStruct_Org.Order
    K_ISE = ISE(K_GMMStruct,GMMStruct_Org);
    fprintf("The ISE is %x \n", K_ISE);
end
