function plot_3D(GMMStruct,x_org,x_cont,title_text,color)

g = figure;

dims = [1,2,3];

plot3(x_org(:,dims(1)),x_org(:,dims(2)),x_org(:,dims(3)),'Marker', '+','LineStyle', 'none','Color', color,'MarkerSize',6);
hold on;
if not(isempty(x_cont))
    plot3(x_cont(:,dims(1)),x_cont(:,dims(2)),x_cont(:,dims(3)),'Marker', 'x','LineStyle', 'none','Color', 'r','MarkerSize',6);
end

FS = 16;

xlabel(sprintf('$x_%d$',dims(1)),'FontSize',FS,'Interpreter','latex');
ylabel(sprintf('$x_%d$',dims(2)),'FontSize',FS,'Interpreter','latex');
zlabel(sprintf('$x_%d$',dims(3)),'FontSize',FS,'Interpreter','latex');

M = GMMStruct.Means; % mean vectors.
R = GMMStruct.Covars; % covariance matrices.

for m = 1:GMMStruct.Order
    plot_gaussian_ellipsoid(M(dims,m).', R(dims,dims,m),2,'cyan');
end

axis tight
grid on

title(title_text);