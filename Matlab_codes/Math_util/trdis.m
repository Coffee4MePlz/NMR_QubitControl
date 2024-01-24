function [sma] = trdis(r1,r2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

m1 = r1-r2;
mm = m1'*m1;
[v,av] = eig(mm);
sm = v*sqrt(real(av))*v';

sma = 0.5*trace(sm);

end

