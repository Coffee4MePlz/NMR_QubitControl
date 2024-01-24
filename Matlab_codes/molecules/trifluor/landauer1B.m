
initspin([1 2 3]);

appulse('F190B',90);

grad(1238,0.001);


appulse('F390B',90);

apdelay(2,3,-70.5);
 
aprefocus('F1180B',0,[1 3]);
aprz(70.5/2,3);
appulse('F390B',270);


grad(1238,0.002);

appulse('F290B',90);

apdelay(2,3,-90);
 
aprefocus('F1180B',0,[1 2]);
aprz(45,2);
appulse('F290B',270);

grad(1238,0.0015);

appulse('F3aB',90);

grad(1238,0.0023);
%---------------------------------------------------------
 appulse('F290B',90);
 
apdelay(2,3,-90);
 
% %---------swap----------------------------------------
%  
  aprefocus('F1180B',0,[1 2]);
%  
 appulse('F1390aB',[90 90]);
% appulse('F190',[90]);
% appulse('F390',[90]);
  apdelay(1,3,90);
%  
  aprefocus('F2180B',0,[1 2]);
%   
  appulse('F1390bB',[270 270]);
% appulse('F190',[270]);
% appulse('F390',[270]);
   apdelay(1,3,90);
   aprefocus('F2180B',0,[1 2]);
%   
  appulse('F1390cB',[0 0]);
% appulse('F190',[0]);
% appulse('F390',[0]);
   apdelay(1,3,90);
   aprefocus('F2180B',0,[1 2]);
% %  
  appulse('F1390dB',[180 180]);
% appulse('F190',[180]);
% appulse('F390',[180]);
% %  
% %  
% %  %------------------------------------
% %  %resetspin([1 3])
% %  
   apdelay(2,3,90);
% % 
   aprefocus('F1180B',0,[1 2]);
% % 
