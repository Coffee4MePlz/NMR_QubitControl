
initspin([1 2 3]);

appulse('F290',90);

grad(1238,0.001);

 
 appulse('F390',270);

apdelay(1,3,70.5);
 
aprefocus('F2180',0,[2 3]);
aprz(-70.5/2,3);
appulse('F390',90);

 
  grad(1238,0.001);

 appulse('F190',270);
 
 apdelay(1,3,90);
  
 aprefocus('F2180',0,[1 2]);
 aprz(-45,1);
 appulse('F190',90);
 
 grad(1238,0.002);
 
appulse('F3a',90);
 
grad(1238,0.0023);
%---------------------------------------------------------
 %appulse('F190',90);
%  
% apdelay(1,3,90);
% aprefocus('F2180',0,[1 2]);
% 
% %---------swap----------------------------------------

% 
%   appulse('F290g',90); appulse('F390g',90);
%   apdelay(2,3,-180);
%   aprefocus('F1180',[0],[1 2]);
%  
%    aprz(-90,2); aprz(-90,3);
%  
%   appulse('F290g',90); appulse('F390g',90);
%   apdelay(2,3,-180);
%   aprefocus('F1180',[0],[1 2]);
% 
%   appulse('F290g',180); appulse('F390g',180);
%   aprz(180,2); aprz(180,3);
%   
%   apdelay(2,3,-180);
%   aprefocus('F1180',[0],[1 2]);
%  
%  
% % %  %------------------------------------

 
 %-----------------------------------------------------------------------
 
 % apdelay(1,3,-90); 
  %aprefocus('F2180',0,[1 2]);
%  
 
global nr;
if nr == 1;  end
if nr == 2; appulse('F290g',0); end
if nr == 3; appulse('F290g',90); end
if nr == 4; appulse('F390g',0);  end
if nr == 5; appulse('F390g',90);  end
if nr == 6; appulse('F290g',0); appulse('F390g',0);  end
if nr == 7; appulse('F290g',90); appulse('F390g',90);  end
if nr == 8; appulse('F290g',0); appulse('F390g',90);  end
if nr == 9; appulse('F290g',90); appulse('F390g',0);  end

appulse('F190',90);


