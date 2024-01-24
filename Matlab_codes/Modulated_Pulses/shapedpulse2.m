function [U roh] = shapedpulse2(shape,pw,phase,ramp,power,roha)
% simulate shapedpulse (shaped pulse on obeserver)
%[U roh] = shapedpulse2(shape,pw,phase,ramp,power,roha);
% pw - > pulse width in microseconds
% shape -> shape of the pulse (like in .RF file) >>>>>> as a matrix <<<<<<<<<
% phase - > phase of the pulse in degree
% ramp ->  introduce a linear ramp on the phase to shift 
%the excitation region (given in Hz)
% power -> it is the maximum amplitude of the pulse (B1 in hz) 
% roha -> density matrix before the pulse
% roh -> density matrix after the pulse
% U ->  Pulse propagator (this incluies the free evoluiton induced by the
% rofs times)


%%%%%%%%%%% obs: whatchout if phase and amp are in reverse order %%%%%%%%%%%

global mol spectro;

phase = phase*pi/180;

shp = shape;

np = length(shp(:,3));  np2 = sum(shp(:,3));

dt = 1e-6*pw/np2;   t = ((1:np)*dt).*shp(:,3)';

wr= 2*pi*power*shp(:,1)/max(shp(:,1));
dphidt = 2*pi*ramp;

ph = shp(:,2)*pi/180 +  phase + dphidt*t' ;
U = 1;

for k=1:np 

         Hrf =   wr(k)*(cos(ph(k))*mol.Ich{1,1} + sin(ph(k))*mol.Ich{1,2});

         U = expm(-1i*(Hrf +  mol.Hint + mol.Hzee)*dt)*U;

end
[Ua roh]= delay(spectro.rof1*1e-6,roha); 
roh = U*roh*U'; 
[Ub roh]= delay(spectro.rof2*1e-6,roh);

U = Ub*U*Ua;
