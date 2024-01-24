
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
aprz(90,3);
aprz(-90,1);

appulse('F190D',0);

apdelay(1,3,180);
aprefocus('F2180D',0,[1 2]);

appulse('F190D',270);
resetspin([1 2 3]);

appulse('F390D',270);
apdelay(1,3,-180);
aprefocus('F2180D',0,[2 3]);

aprz(90,3); %theta/2
appulse('F390D',90);
aprz(90,1); %theta/2
aprz(90,3);
aprz(-90,1);

appulse('F190D',0);

apdelay(1,3,180);
aprefocus('F2180D',0,[1 2]);
appulse('F190D',270);


% 
%  appulse('F1390aD',[90 90]);
%  apdelay(1,3,180);
%   
%   aprefocus('F2180D',0,[1 2]);
%    
%   appulse('F1390bD',[270 270]);
% 
%   apdelay(1,3,180);
%   aprefocus('F2180D',0,[1 2]);
%    
%   appulse('F1390cD',[0 0]);
%   apdelay(1,3,180);
%   aprefocus('F2180D',0,[1 2]);
% %   
 appulse('F1390dD',[0 0]);
 %appulse('F1390dD',[180 180]);
 
 
 