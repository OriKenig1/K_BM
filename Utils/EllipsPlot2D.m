function EllipsPlot2D(m,C,binnum,c)



% This function plots a 2-dimensional ellipsoid corresponding to a 2-dimensional Gaussian

% random vector.



% Usage:

% EllipsPlot2D(m,C,binnum,c)



% Inputs:

% m - 2x1 mean vector.

% C - 2x2 covariance matrix.

% binnum - number of plotting points.

% c - color of plot (string).



% Version 1.0, Mai 2006.

% By Koby Todros  todros@widemed.com





[V,D] = eig(C);

r = 2*sqrt(diag(D));

temp = m;

m = [0;0];

x = m(1)+linspace(-r(1),r(1),binnum);

y1 = real(m(2)+r(2)*sqrt(1-((x-m(1))/r(1)).^2));

y2 = real(m(2)-r(2)*sqrt(1-((x-m(1))/r(1)).^2));



m = temp;

P = [x',y1'];

P = P*V;

P(:,1) = P(:,1)+m(1);

P(:,2) = P(:,2)+m(2);

plot(P(:,1),P(:,2),c,'LineWidth',1);

hold on

P = [x',y2'];

P = P*V;

P(:,1) = P(:,1)+m(1);

P(:,2) = P(:,2)+m(2);

plot(P(:,1),P(:,2),c,'LineWidth',1);

hold on 
plot(m(1),m(2),'x');


