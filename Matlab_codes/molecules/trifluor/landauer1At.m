
initspin([1 2 3]);

appulse('F190A',90);

grad(1238,0.001);


appulse('F390A',90);

apdelay(2,3,-70.5);
 
aprefocus('F1180A',0,[1 3]);
aprz(70.5/2,3);
appulse('F390A',270);


grad(1238,0.002);

appulse('F290A',90);

apdelay(2,3,-90);
 
aprefocus('F1180A',0,[1 2]);
aprz(45,2);
appulse('F290A',270);

grad(1238,0.0015);

appulse('F3aA',90);

grad(1238,0.0023);
%---------------------------------------------------------
aprz(90,3);
aprz(-90,1);

appulse('F190A',0);

apdelay(1,3,180);
aprefocus('F2180A',0,[1 2]);

appulse('F190A',270);
resetspin([1 2 3]);

appulse('F390A',270);
apdelay(1,3,-180);
aprefocus('F2180A',0,[2 3]);

aprz(90,3); %theta/2
appulse('F390A',90);
aprz(90,1); %theta/2
aprz(90,3);
aprz(-90,1);

appulse('F190A',0);


apdelay(1,3,180);
aprefocus('F2180A',0,[1 2]);
appulse('F190A',270);


%   appulse('F1390aA',[90 90]);
%  
%   apdelay(1,3,180);
%    
%   aprefocus('F2180A',0,[1 2]);
%      
%    appulse('F1390bA',[270 270]);
%   
%    apdelay(1,3,180); 
%    aprefocus('F2180A',0,[1 2]);
%      
%    appulse('F1390cA',[0 0]);
%     
%    apdelay(1,3,180);
%     
%    aprefocus('F2180A',0,[2 3]);
%     
%    appulse('F1390dA',[180 180]);
    