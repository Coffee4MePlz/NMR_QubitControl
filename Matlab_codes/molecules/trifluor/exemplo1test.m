
%consider the spins on Z direction
resetspin([1 2 3]); 

%ask the pi/2 pulse using shape F290
appulse('F290',0);  

%ask the delay 
apdelay(1,2,180);

%ask a refocising pulse to cancel the J23 
aprefocus('F3180',0,[2 3]);

