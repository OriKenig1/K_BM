clear 
close all

Dimension = 10;
Order = 3;
N = 300;

MaxIter = 50;
tol = 10^-13;

dim_vec = 1:Dimension;

% -------------------------------------------
% Original PDF

v1 = [1,-1,1,-1,1,-1,1,-1,1,-1]';

%mu_1 = [10,-2,10*ones(1,Dimension-2)]';
mu_1 = 2*ones(1,Dimension)';
%mu_2 = 3*ones(Dimension,1);
mu_2 = 14*ones(1,Dimension)';
% mu_3 = 13*ones(Dimension,1);
mu_3 = mu_1 + v1*8;
mu_Org = [mu_1,mu_2,mu_3];

Sigma_1 = 5 * 0.95 .^ (abs( dim_vec'-dim_vec)) ;
Sigma_2 = 3 * eye(Dimension);
Sigma_3 = 3 * 0.8 .^ (abs( dim_vec'-dim_vec)) ;
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
mu_2 = -40*ones(1,Dimension)' + v1*5;
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

epsilon = 0.02;

% -------------------------------------------
% Initialize 

[x,x_org,x_cont] = Generate_Obs(GMMStruct_Org,GMMStruct_cont,N,epsilon);

GMMStruct_0 = KMEANS_initGuess(x,Order);

I = 1:0.5:20;

% ----------------------------------------- %
% Plot Original

fh = figure();
fh.WindowState = 'maximized';

for j = 1:Dimension/2

d1 = 2*(j-1) + 1;
d2 = 2*j;

subplot(3,2,j);
plot(x_org(:,d1),x_org(:,d2),'b+','MarkerSize',2);
hold on;
plot(x_cont(:,d1),x_cont(:,d2),'r*','MarkerSize',2);
xlabel('$x_1$','FontSize',16,'Interpreter','latex');
ylabel('$x_2$','FontSize',16,'Interpreter','latex');
for i = 1:Order
    hold on
    EllipsPlot2D((GMMStruct_Org.Means(d1:d2,i)),(GMMStruct_Org.Covars(d1:d2,d1:d2,i)),100,'black');
    plot((GMMStruct_Org.Means(d1,i)),(GMMStruct_Org.Means(d2,i)),'blackx','MarkerSize',2);
  
end

title(sprintf("Dimensions %d - %d",d1,d2));

end
sgtitle("Original");
% ----------------------------------------- %

% ----------------------------------------- %
% Run K-BM %

[GMMStruct,h_opt] = K_BM(x,Order,I,GMMStruct_0);
K_MISE = MISE(GMMStruct_Org,GMMStruct);
display(K_MISE,'K-BM MISE');
display(h_opt,'h_opt');

% ----------------------------------------- %
% Plot K-BM

fh = figure();
fh.WindowState = 'maximized';

for j = 1:Dimension/2

d1 = 2*(j-1) + 1;
d2 = 2*j;

subplot(3,2,j);
p1 = plot(x_org(:,d1),x_org(:,d2),'b+','MarkerSize',2);
hold on;
p2 = plot(x_cont(:,d1),x_cont(:,d2),'r*','MarkerSize',2);

xlabel('$x_1$','FontSize',16,'Interpreter','latex');
ylabel('$x_2$','FontSize',16,'Interpreter','latex');
for i = 1:Order
    hold on
    EllipsPlot2D((GMMStruct.Means(d1:d2,i)),(GMMStruct.Covars(d1:d2,d1:d2,i)),100,'black');
    plot((GMMStruct.Means(d1,i)),(GMMStruct.Means(d2,i)),'blackx','MarkerSize',2);
  
end

title(sprintf("Dimensions %d - %d",d1,d2));

end
sgtitle("K-BM");
% ----------------------------------------- %


