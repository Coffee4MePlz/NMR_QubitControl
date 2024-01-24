
%consider the spins on Z direction
 initspin([1 2 3])
 
 appulse('F190',90);  
 
 aprz(-49.1,1);
 
 apdelay(1,2,98.2);
 
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
  
  aprz(67.8,3);
  
  apdelay(2,3,-135.6);
  
  aprefocus('F1180',0,[2 3]);
  
  appulse('F390',90); 
  
   grad(1234,0.0014);
 
  appulse('F290',270);  
  
  aprz(45,2);
  
  apdelay(2,3,-90);
  
  aprefocus('F1180',0,[1 2]);
  
  appulse('F290',90); 
  
  grad(1234,0.0023);
appulse('F290',90); 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% unresetspin([1 3]);
 
%appulse('F290',90);
%appulse('F390',270);
% 
% appulse('F2390',[90 270]);
% 
%  apdelay(1,3,-180);
%  
%  aprefocus('F2180',0,[2 3]);
% 
%  appulse('F390',0); 
%  aprz(90,1);
%  aprz(90,3);
%  
%  apdelay(1,2,-180);
%  
%  aprefocus('F3180',0,[1 3]);
%  
%  appulse('F190',180);  
%  
%  apdelay(1,3,-90);
%  
%  aprefocus('F2180',0,[1 2]);
%  
%  appulse('F190',0);  
%  
%  apdelay(1,2,180);
%  aprefocus('F3180',0,[1 3]);
%   appulse('F190',90);  
%  
%  apdelay(1,2,-90);
%  opthere;
%  
%  apdelay(1,3,-90);
%  
%  aprefocus('F2180',0,[1 2]);
%  
%  aprz(45,1);
%  appulse('F190',270);  
%  
%  apdelay(2,3,90);
%  aprefocus('F1180',0,[1 3]);
%  
%  aprz(-45,2);
%  aprz(-45,3);
%  
%  appulse('F390',270);  
%  apdelay(1,3,-180);
%  
%  aprefocus('F2180',0,[1 2]);
%  
%  appulse('F390',0); 
%  aprz(90,1);
%  aprz(90,3);
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %appulse('F290',90); 
%  
%  
%  
%  
%  
 
 
 
 
 
 
