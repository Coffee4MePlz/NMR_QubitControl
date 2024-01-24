%%%%%%%% CUSTOMCOLORMAP ------------------------------------------------------------------------------------------------------------------------
%%%%%%%% function [myColorMap] = custoncolormap(flag,nfig1);
%%%%%%%% Code to create a custom 28 level colormap, then rescale it to 256 levels to be smoother.
%%%%%%%% TO MODIFY THE MAP: Choose a diferent set of exadecimal colors in
%%%%%%%% hexMap array
%%%%%%%% flag 0 nofigure | flag = 1 colormap figure | nfig1 number of the
%%%%%%%% figure for the color map
function [myColorMap] = custoncolormap(flag,nfig1);
format long g;
format compact;
fontSize = 18; % Fontsize for the colormap figure
%fprintf('Beginning to run %s.m ...\n', mfilename);
%echo off;
hexMap = {'C0C0C0', '808080', '404040', '000000', 'FF99CC', '9999FF', '3333FF', '000099', '3399FF', '0066CC', '99CCFF', '66B2FF', '66FFFF', '006633', '00CC66', '66FF66', '00FF00', '009900', 'FFFF99', 'FFFF00', 'CCCC00', 'FFB266', 'CC6600', '994C00', 'FF9999', 'FF0000', 'CC0000', '990000'};
myColorMap = zeros(length(hexMap), 3); % Preallocate varable
numLevels = length(hexMap);
r = zeros(numLevels, 1);
g = zeros(numLevels, 1);
b = zeros(numLevels, 1);
for k = 1 : numLevels
	thisCell = hexMap{k};
	r(k) = hex2dec(thisCell(1:2));
	g(k) = hex2dec(thisCell(3:4));
	b(k) = hex2dec(thisCell(5:6));
end
% Now that the individual color channel mappings have been created,
% stitch the red, green, and blue vectors together to create the 256-by-3 colormap matrix.
myColorMap = [r, g, b];
myColorMap = myColorMap / 255; % Normalize to range 0-1
% Scale to 256 gray levels, instead of the current 28.
numLevels = 256;
myColorMap = imresize(myColorMap, [numLevels, 3]);
% Sometimes it exceeds 1 so rescale 0-1
myColorMap = rescale(myColorMap, 0, 1);
%===============================================================================
% Plot the red, green, and blue components of the custom colormap so we can see what it looks like.
if flag == 1
    figure(nfig1)
    r = myColorMap(:, 1);
    g = myColorMap(:, 2);
    b = myColorMap(:, 3);
    igs = 0 : length(r) - 1; % Input gray scale axis.
    plot(igs, r, 'r-', 'LineWidth', 3);
    hold on;
    plot(igs, g, 'g-', 'LineWidth', 3);
    plot(igs, b, 'b-', 'LineWidth', 3);
    grid on;
    caption = sprintf('Color Map with %d Rows', length(r));
    title(caption, 'FontSize', fontSize);
    xlabel('Input Gray Level', 'FontSize', fontSize);
    ylabel('Output Color Channel Level', 'FontSize', fontSize);
    xlim([min(igs), max(igs)]);
end

end
%===============================================================================
