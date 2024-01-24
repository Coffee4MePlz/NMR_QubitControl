%%%%%%%% PrettyPlotMatrix 2 ------------------------------------------------------------------------------------------------------------------------
%%%%%%%% Pretty Plot of the density matrix with the argument as colour map in one figure 
%%%%%%%% Imputy aaroh = density matrix | nunfig = figure number  
function PrettyPlotMatrix2(aaroh,nunfig);

transp = 0.6; % opacidade das colunas
fontsizetitle = 19; % tamanho da fonte titulo
fontsizeaxes = 17; % tamanho da fonte eixos

[myColorMap] = custoncolormap(1,nunfig+3); % gera um color map customizado
[myColorMap] = colormap(jet);

figure(nunfig)

Y = aaroh;
h = bar3(abs(real(Y))); % graficos de barras para a parte real

% Construindo o colormap para o argumento de rho

% colormap jet % escolhe o mapa de cor do matlab a ser usado tambem pode ser personalizado
% cm = get(gcf,'colormap');  % Atribui o colormap do matlab ao array cm

cm = myColorMap;  % Atribui o colormap personalizado ao array cm

for ii = 1:length(Y); 
    for jj = 1:length(Y);
        Yarg_color = 1+floor(255*(angle(Y(ii,jj)) + pi)/(2*pi)); % Retorna argumeto compexo de Y(ii,jj) que varia de 1 (-pi) a 256 (pi) 
        cm2{ii,jj} = cm(Yarg_color,:);
    end; 
end;

% Recoloring the bar plot
cnt = 0;
for jj = 1:length(h)
    xd = get(h(jj),'xdata');
    yd = get(h(jj),'ydata');
    zd = get(h(jj),'zdata');
    delete(h(jj))    
    idx = [0;find(all(isnan(xd),2))];
    for ii = 1:length(idx)-1
         cnt = cnt + 1;
         S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
                         yd(idx(ii)+1:idx(ii+1)-1,:),...
                         zd(idx(ii)+1:idx(ii+1)-1,:),...
                         'facecolor',cm2{ii,jj},'FaceAlpha',transp);
     end
end

colormap(myColorMap); % define o colormap para a barra de cores
colorbar;      % Plota a barra de corres

colorbar('Ticks',[0,0.125,0.25,0.375,0.5],...
         'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})

zlim([-0.05 0.6])
yticklabels({'00','01','10','11'})
xticklabels({'00','01','10','11'})
ax = gca(figure(nunfig));
ax.FontSize = fontsizeaxes;

% t=title(['Arg(\rho) as colours']);
% t.FontSize = fontsizetitle;

end