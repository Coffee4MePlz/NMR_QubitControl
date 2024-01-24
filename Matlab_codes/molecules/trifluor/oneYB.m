
%consider the spins on Z direction
 initspin([1 2 3])
 
 appulse('F190',90);  
 

 %aprz(-49.1,1);
 aprz(-49.106605350869096,1);
 
 %apdelay(1,2,98.2);
 apdelay(1,2,98.213210701738191);
 
 
 aprefocus('F3180',0,[1 3]);
 
 
 appulse('F190',270); 
 
 grad(1234,0.001);
 
  appulse('F390',270);  
  
  aprz(90,3);
  
  apdelay(2,3,-180);
  
  aprefocus('F1180',0,[1 3]);
  
  appulse('F390',90); 
  
   grad(1234,0.0017);
  
  appulse('F290',270);  
  
  aprz(45,2);
  
  apdelay(1,2,-90);
  
  aprefocus('F3180',0,[2 3]);
  
  appulse('F290',90); 
  
  grad(1234,0.0022);
 
  appulse('F390',270);  
  

  %aprz(67.8,3);
  aprz(67.792345701403519,3);
  
  %apdelay(2,3,-135.6);
  apdelay(2,3,-135.5846914028070);
  
  
  aprefocus('F1180',0,[2 3]);
  
  appulse('F390',90); 
  
   grad(1234,0.0014);
 
  appulse('F290',270);  
  
  aprz(45,2);
  
  apdelay(2,3,-90);
  
  aprefocus('F1180',0,[1 2]);
  
  appulse('F290',90); 
  
  grad(1234,0.0023);

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %==========================================================================
% 
%   
% %=================================YB=======================================

%      appulse('F190g',10);
%      appulse('F190g',190);
%      appulse('F190g',10);
%      appulse('F390g',90);
%   appulse('F390g',90);
%   appulse('F190g',90);
%   appulse('F390g',270);
%   appulse('F190g',270);
%   appulse('F390g',90);
%   appulse('F190g',90);
%  
%  
   appulse('F190g',270);
   appulse('F190g',90);
   appulse('F190g',270);
% % %  
% % %  
   
   appulse('F390g',90);
   appulse('F390g',270);
   appulse('F390g',90);
% %   
%  
%  aprz(90,1);
%  appulse('F190g',180);
%  
%  aprz(90,3);
%  appulse('F390g',180);
% 
% 
   
%  
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Had%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
  appulse('F290g',90);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CSWAP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
 
  appulse('F190g',90);
 
  apdelay(1,3,180);
  
  aprefocus('F2180',0,[1 2]);
 
  appulse('F190',0); 
  aprz(-90,1);
  aprz(90,3);
  
  
  %---
  
  apdelay(1,3,-180);
  
  aprefocus('F2180',0,[2 3]);
  
  appulse('F390',180);  
  
  apdelay(2,3,-90);
  
  aprefocus('F1180',0,[1 3]);
  
  appulse('F390',180);  
  
  apdelay(1,3,-180);
  aprefocus('F2180',0,[2 3]);
   appulse('F390',90);  
  aprz(180,3);
  
  apdelay(2,3,-90);
  aprefocus('F1180',0,[1 2]);
  opthere;
  
  apdelay(1,3,-90);
  
  aprefocus('F2180',0,[2 3]);
  
  aprz(45,3);
  appulse('F390',270);  
  
  apdelay(1,2,90);
  aprefocus('F3180',0,[1 3]);
  
  aprz(-45,1);
  aprz(-45,2);
  

  appulse('F190',270);

  apdelay(1,3,-180);
  
  aprefocus('F2180',0,[1 2]);
 
  appulse('F190g',0);
  aprz(-180,1);
  aprz(-90,1);
  aprz(90,3);
  
  
  
  
%  appulse('F290g',90);
% 
% 
% 
%   appulse('F190g',90);
%  
%   apdelay(1,3,180);
%   
%   aprefocus('F2180',0,[1 2]);
%  
%   appulse('F190',180); 
%   aprz(-90,1);
%   aprz(90,3);
%   
%   apdelay(1,3,180);
%   
%   aprefocus('F2180',0,[2 3]);
%   
%   appulse('F390',180);  
%   
%   apdelay(2,3,-90);
%   
%   aprefocus('F1180',0,[1 3]);
%   
%   appulse('F390',0);  
%   
%   apdelay(1,3,-180);
%   aprefocus('F2180',0,[2 3]);
%    appulse('F390',90);  
%   %aprz(180,3);
%   
%   apdelay(2,3,-90);
%   aprefocus('F1180',0,[1 2]);
%   opthere;
%   
%   apdelay(1,3,90);
%   
%   aprefocus('F2180',0,[2 3]);
%   
%   aprz(45,3);
%   appulse('F390',270);  
%   
%   apdelay(1,2,-90);
%   aprefocus('F3180',0,[1 3]);
%   
%   aprz(-45,1);
%   aprz(-45,2);
%   aprz(-270,1);
%   
%   appulse('F190',90);
%   
%   
%   apdelay(1,3,-180);
%   
%   aprefocus('F2180',0,[1 2]);
%  
%   appulse('F190g',0);
%   aprz(-180,1);
%   aprz(-90,1);
%   aprz(90,3);



















  
   % appulse('F290g',90);
    
    