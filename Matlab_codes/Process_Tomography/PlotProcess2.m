%%%%%%%% PlotProcess ------------------------------------------------------------------------------------------------------------------------
%%%%%%%% Pretty Plot of real and imaginary part the the density matrix in two separate figures 
function PlotProcess(Choimat,nunfig);

%%%%%%% Representacao do roh
f = figure(nunfig);
f.Position = [1200 600 720 300];
%subplot(1,2,1)

transp = 0.6; % opacidade das colunas

PP = real(vpa(Choimat));
b = bar3(PP);
%colorbar;

for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% For each surface in bar graph - Loop through
% Face Alpha controls transparency of surface
% 0 is completely transparent while 1 is opaque
for i = 1:length(b)
  b(i).FaceAlpha = transp;
end

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis

yticklabels({'$i 1$','$\sigma_x$','$\sigma_y$','$\sigma_z$'});
xticklabels({'$i 1$','$\sigma_x$','$\sigma_y$','$\sigma_z$'});
title('Re')

% xlabel('i', 'FontSize', 14);
% ylabel('j', 'FontSize', 14);
% zlabel('p_{i|j}', 'FontSize', 14);

ax = gca(figure(nunfig));
ax.FontSize = 14;

%{
% Imaginary part
PP = imag(Choimat);
teste =  all(all(PP == 0) == 1);

subplot(1,2,2)

b = bar3(PP);
%colorbar;

for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% For each surface in bar graph - Loop through
% Face Alpha controls transparency of surface
% 0 is completely transparent while 1 is opaque
for i = 1:length(b)
  b(i).FaceAlpha = transp;
end

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex lfor y-axis

zlim([-1.005 1.005]);
yticklabels({'$i 1$','$\sigma_x$','$\sigma_y$','$\sigma_z$'});
xticklabels({'$i 1$','$\sigma_x$','$\sigma_y$','$\sigma_z$'});

ax = gca(figure(nunfig));
ax.FontSize = 14;
title('Im')

% Valores sobre as barras
% [X,Y] = meshgrid(1:size(PP,2), 1:size(PP,1));
% text(X(:), Y(:), PP(:), num2str(PP(:)), 'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize', 13)

% xlabel('i', 'FontSize', 14);
% ylabel('j', 'FontSize', 14);
% zlabel('p_{i|j}', 'FontSize', 14);
%}
end