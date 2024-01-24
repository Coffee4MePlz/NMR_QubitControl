%%%%%%%% PrettyPlotMatrix ------------------------------------------------------------------------------------------------------------------------
%%%%%%%% Pretty Plot of real and imaginary part the the density matrix in two separate figures 
function PrettyPlotMatrix(aaroh,nunfig);

%%%%%%% Representacao do roh
figure(nunfig)
transp = 0.6; % opacidade das colunas

PP = real(aaroh);
b = bar3(PP);
colorbar;

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

yticklabels({'00','01','10','11'});
xticklabels({'00','01','10','11'});

% Valores sobre as barras
% [X,Y] = meshgrid(1:size(PP,2), 1:size(PP,1));
% text(X(:), Y(:), PP(:), num2str(PP(:)), 'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize', 13)

% xlabel('i', 'FontSize', 14);
% ylabel('j', 'FontSize', 14);
% zlabel('p_{i|j}', 'FontSize', 14);

ax = gca(figure(nunfig));
ax.FontSize = 14;

% Imaginary part
figure(nunfig+1);
PP = imag(aaroh);
b = bar3(PP);
colorbar;

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

zlim([-0.505 0.505]);
yticklabels({'00','01','10','11'});
xticklabels({'00','01','10','11'});

ax = gca(figure(nunfig+1));
ax.FontSize = 14;

% Valores sobre as barras
% [X,Y] = meshgrid(1:size(PP,2), 1:size(PP,1));
% text(X(:), Y(:), PP(:), num2str(PP(:)), 'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize', 13)

% xlabel('i', 'FontSize', 14);
% ylabel('j', 'FontSize', 14);
% zlabel('p_{i|j}', 'FontSize', 14);

end
