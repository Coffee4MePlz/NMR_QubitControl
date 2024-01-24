function [r th power np dt]=composepul(shps,pw,ang,phase,targ,typs);
  % usage --> [wr{k} ph{k} aa np dt] = composepul(shp,pw,ang2,phase2,targ2,typs2);

% shps := {{'pulse type.RF'}} cell array % in the format of .RF [phase,amplitude, timestep]
% pw := float % pulse width,     %  targ := vector of int % target qubits
% ang := float / maybe an array (?) % angle of rotation
% typs := cell array {{'n'}} for normal (?) or {{'g'}} for gradient (?)
%r and th are polar coordinates

global mol;

rr = 0; ii = 0;

for k=1:length(shps)
    shp = shps{k};
    np = length(shp(:,3));  np2 = sum(shp(:,3));
    dt = pw/np2;   t = ((1:np)*dt).*shp(:,3)';
    
    if typs{k} == 'n'
    
        ramp(k) = 2*pi*mol.dq(targ(k)); %ramp is zero where the pulse is focused
        pp = pi*(shp(:,1) - shp(1,1))/180;
        wr = pi*ang(k)*abs(shp(:,2)/sum(cos(pp).*shp(:,2)*dt*180)); % (?) whats the meaning of this ? and if cos(pp)=0?
        
    else
        ramp(k) = 0;         
        wr = shp(:,2);
    
    end
    
    ph = shp(:,1)*pi/180 + phase(k)*pi/180 + ramp(k)*t';
    
    rr = rr + wr.*cos(ph);
    ii = ii + wr.*sin(ph);

end

[th r] = cart2pol(rr,ii); %converts cartesian to polar coordinates

power = max(r)/(2*pi);
