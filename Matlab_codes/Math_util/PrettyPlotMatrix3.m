%%%%%%%% PrettyPlotMatrix ------------------------------------------------------------------------------------------------------------------------
%%%%%%%% Pretty Plot of real and imaginary part the the density matrix in two separate figures 

% now working for N Qbits

function PrettyPlotMatrix3(aaroh,nunfig,flag);

% defining size of rho
[p1,q1] = size(aaroh);

%%%%%%% Representacao do roh
f = figure(nunfig);
f.Position = [1200 600 1040 400];
tiledlayout(1,2)
nexttile;


transp = 0.6; % opacidade das colunas

PP = real(aaroh);
b = bar3(PP);
colorbar;
%clim([-6, 6]);
zlim([-1.1 1.1])
clim([-1, 1]);
%zlim([0 1]);

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

set_ticklabels(p1);

% Valores sobre as barras
if flag == 1
    [X,Y] = meshgrid(1:size(PP,2), 1:size(PP,1));
    %aux = round(PP,4);
    text(X(:), Y(:), PP(:), num2str(round(PP(:),4)), 'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize', 9)
end

% xlabel('i', 'FontSize', 14);
% ylabel('j', 'FontSize', 14);
% zlabel('p_{i|j}', 'FontSize', 14);

ax = gca(figure(nunfig));
%ax.FontSize = 14;
set_FontSize(p1,ax);
title("Real(rho)");
% Imaginary part
nexttile;
%figure(nunfig+1);


PP = imag(aaroh);
b = bar3(PP);
colorbar;
%clim([-5, 5]);
clim([-0.5, 0.5]);
%zlim([-5 5]);
zlim([-0.505 0.505])

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

%ax = gca(figure(nunfig+1));
ax = gca(figure(nunfig));
%ax.FontSize = 14;
set_ticklabels(p1)
set_FontSize(p1,ax);
title("Im(rho)");
% Valores sobre as barras
if flag == 1
    [X,Y] = meshgrid(1:size(PP,2), 1:size(PP,1));
    text(X(:), Y(:), PP(:), num2str(round(PP(:),4)), 'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize', 9)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set_ticklabels(size)
    BitStr = {};
    ticknum = ones(1,size);
    str_size = log2(size);
    for j= 1:size
        BitStr{j} = dec2bin(j-1,str_size);
        ticknum(j) = j;
    end
    yticks(ticknum);
    xticks(ticknum);
    xticklabels(BitStr);
    yticklabels(BitStr);
end

function set_FontSize(size,ax)
    if size <8
        ax.FontSize = 14;
    elseif size == 8
        ax.FontSize = 11;
    elseif size > 8
        ax.FontSize = 6;
    end

end
