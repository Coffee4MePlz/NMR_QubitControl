
initspin([1 2 3]);

appulse('F190D',90);

grad(1238,0.001);


appulse('F390D',90);

apdelay(2,3,-70.5);
 
aprefocus('F1180D',0,[1 3]);
aprz(70.5/2,3);
appulse('F390D',270);


grad(1238,0.002);

appulse('F290D',90);

apdelay(2,3,-90);
 
aprefocus('F1180D',0,[1 2]);
aprz(45,2);
appulse('F290D',270);

grad(1238,0.0015);

appulse('F3aD',90);

grad(1238,0.0023);
%---------------------------------------------------------
 appulse('F290D',90);
 
apdelay(2,3,-90);
 
% %---------swap----------------------------------------
%  
  aprefocus('F1180D',0,[1 2]);
    
  appulse('F2180D',90);
  
  apdelay(2,3,-90);
% % 
   aprefocus('F1180D',0,[1 2]);
% % 
