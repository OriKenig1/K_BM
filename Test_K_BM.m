clear 
close all

% -------------------------------------------
% Tune this parameters to your liking 

Dimension = 10;
Order = 3;
N = 300;

MaxIter = 50;
tol = 10^-13;

epsilon = 0.1;

I = 1:0.5:20;

d1 = 1;
d2 = 2;

v1 = [1,-1,1,-1,1,-1,1,-1,1,-1]';

% -------------------------------------------

% -------------------------------------------
% Original PDF

dim_vec = 1:Dimension;

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

alpha_Org = [0.3,0.3,0.4];

GMMStruct_Org = CreateGMMStruct(Dimension, Order);
GMMStruct_Org.Alpha = alpha_Org;
GMMStruct_Org.Means = mu_Org;
GMMStruct_Org.Covars = Sigma_Org;

% -------------------------------------------

% -------------------------------------------
% Contamination PDF

Order_cont = 2;

mu_1 = 10*ones(1,Dimension)' - v1*10;
mu_2 = -20*ones(1,Dimension)' + v1*5;
mu_cont = [mu_1,mu_2];

Sigma_1 = 10 * 0.9 .^ (abs( dim_vec'-dim_vec));
Sigma_2 = 10 * 0.9 .^ (abs( dim_vec'-dim_vec));
Sigma_cont(:,:,1) = Sigma_1;
Sigma_cont(:,:,2) = Sigma_2;

alpha_cont = [0.5,0.5];

GMMStruct_cont = CreateGMMStruct(Dimension, Order_cont);
GMMStruct_cont.Alpha = alpha_cont;
GMMStruct_cont.Means = mu_cont;
GMMStruct_cont.Covars = Sigma_cont;

% -------------------------------------------
% Initialize 

% Use this to generate new samples
%[x,x_org,x_cont] = Generate_Obs(GMMStruct_Org,GMMStruct_cont,N,epsilon);

load('samples.mat')

GMMStruct_0 = KMEANS_initGuess(x,Order);

% ----------------------------------------- %
% Plot Original

figure();

plot(x_org(:,d1),x_org(:,d2),'go','MarkerSize',3);
hold on;
plot(x_cont(:,d1),x_cont(:,d2),'r*','MarkerSize',10);
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');

M = [ GMMStruct_Org.Means(d1,:); GMMStruct_Org.Means(d2,:) ];
for i = 1:Order
    R = [ GMMStruct_Org.Covars(d1,d1,i), GMMStruct_Org.Covars(d1,d2,i); ...
        GMMStruct_Org.Covars(d2,d1,i), GMMStruct_Org.Covars(d2,d2,i)];
    hold on
    EllipsPlot2D(M(:,i),R,100,'black');
    plot(M(1,i),M(2,i),'blackx','MarkerSize',10);
  
end

title(sprintf("Original : Dimensions %d - %d",d1,d2));
ax = gca;
ax.FontSize = 14;

% ----------------------------------------- %

% ----------------------------------------- %
% Run K-BM %

[GMMStruct,h_opt] = K_BM(x,Order,I,GMMStruct_0);
K_ISE = ISE(GMMStruct_Org,GMMStruct);
display(K_ISE,'K-BM ISE');
display(h_opt,'h_opt');

% ----------------------------------------- %
% Plot K-BM

figure();

plot(x_org(:,d1),x_org(:,d2),'go','MarkerSize',3);
hold on;
plot(x_cont(:,d1),x_cont(:,d2),'r*','MarkerSize',10);
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');

M = [ GMMStruct.Means(d1,:); GMMStruct.Means(d2,:) ];
for i = 1:Order
    R = [ GMMStruct.Covars(d1,d1,i), GMMStruct.Covars(d1,d2,i); ...
        GMMStruct.Covars(d2,d1,i), GMMStruct.Covars(d2,d2,i)];
    hold on
    EllipsPlot2D(M(:,i),R,100,'black');
    plot(M(1,i),M(2,i),'blackx','MarkerSize',10);
  
end

title(sprintf("K-BM : Dimensions %d - %d",d1,d2));
ax = gca;
ax.FontSize = 14;

% ----------------------------------------- %