
initspin([1 2 3]);

appulse('F190C',90);

grad(1238,0.001);


appulse('F390C',90);

apdelay(2,3,-70.5);
 
aprefocus('F1180C',0,[1 3]);
aprz(70.5/2,3);
appulse('F390C',270);


grad(1238,0.002);

appulse('F290C',90);

apdelay(2,3,-90);
 
aprefocus('F1180C',0,[1 2]);
aprz(45,2);
appulse('F290C',270);

grad(1238,0.0015);

appulse('F3aC',90);

grad(1238,0.0023);
%---------------------------------------------------------
 appulse('F290C',90);
 
apdelay(2,3,-90);
 
% %---------swap----------------------------------------
%  
  aprefocus('F1180C',0,[1 2]);
%  
 appulse('F1390aC',[90 90]);
% appulse('F190',[90]);
% appulse('F390',[90]);
  apdelay(1,3,30);
%  
  aprefocus('F2180C',0,[1 2]);
%   
  appulse('F1390bC',[270 270]);
% appulse('F190',[270]);
% appulse('F390',[270]);
   apdelay(1,3,30);
   aprefocus('F2180C',0,[1 2]);
%   
  appulse('F1390cC',[0 0]);
% appulse('F190',[0]);
% appulse('F390',[0]);
   apdelay(1,3,30);
   aprefocus('F2180C',0,[1 2]);
% %  
  appulse('F1390dC',[180 180]);
% appulse('F190',[180]);
% appulse('F390',[180]);
% %  
% %  
% %  %------------------------------------
% %  %resetspin([1 3])
% %  
   apdelay(2,3,90);
% % 
   aprefocus('F1180C',0,[1 2]);
% % 
