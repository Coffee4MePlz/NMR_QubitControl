global nr
initspin([1 2 3]);

appulse('F190',90);

grad(1238,0.001);


appulse('F390',90);

apdelay(2,3,-70.5);
 
aprefocus('F1180',0,[1 3]);
aprz(70.5/2,3);
appulse('F390',270);


grad(1238,0.002);

appulse('F290',90);

apdelay(2,3,-90);
 
aprefocus('F1180',0,[1 2]);
aprz(45,2);
appulse('F290',270);

 grad(1238,0.0015);
 
 appulse('F3a',90);
 
 grad(1238,0.0023);
%---------------------------------------------------------


if nr == 1
	appulse('F290',90);
end

if nr == 2
	appulse('F390',0);  
    apdelay(1,3,180);
    aprefocus('F2180',0,[2 3]);
    appulse('F390',270);
    resetspin([3]);
	appulse('F290',90);
end


if nr == 3
	
	appulse('F290',90);
    apdelay(2,3,-720);
    aprefocus('F1180',0,[1 2]);
 
	
	appulse('F390',0);  
    apdelay(1,3,180);
    aprefocus('F2180',0,[2 3]);
    appulse('F390',270);
    resetspin([3]);
	
	apdelay(2,3,720);
    aprefocus('F1180',0,[1 2]);
end

% 
% if nr ==1;  end
% if nr ==2;    appulse('F190',0); end
% if nr ==3;   appulse('F190',90); end
% if nr ==4;   appulse('F190',0); appulse('F390',0); end