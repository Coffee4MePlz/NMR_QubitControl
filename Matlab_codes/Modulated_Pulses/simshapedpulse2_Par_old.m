function [U, roh] = simshapedpulse2_Par(shape_obs,shape_dec,pw_obs,pw_dec,phase_obs,phase_dec,... 
    ramp_obs,ramp_dec,power_obs,power_dec,roha, L_Spec, L_Mol)
% simulate simshapedpulse (shaped pulse simultaneously on obeserver and decopler channels)
% [U, roh] = simshapedpulse(shape_obs,shape_dec,pw_obs,pw_dec,phase_obs,phase_dec,ramp_obs,ramp_dec,power_obs,power_dec,roha)
% pw - > pulse width in microseconds
% shape -> shape of the pulse (like in .RF file without any comment)
% phase - > phase of the pulse in degree
% ramp ->  introduce a linear ramp on the phase to shift the excitation region (given in Hz)
% power -> it is the maximum amplitude of the pulse (B1 in Hz) 
% roha -> density matrix before the pulse
% roh -> density matrix after the pulse
% U ->  Pulse propagator (this incluies the free evoluiton induced by the rofs times)
% Both shaps (obs and dec) must have the same points numbers

phase_obs = phase_obs*pi/180; % initial phase in radians
phase_dec = phase_dec*pi/180; % initial phase in radians

shp_obs = shape_obs;
shp_dec = shape_dec;

np = length(shp_obs(:,3));  np2 = sum(shp_obs(:,3));

dt = 1e-6*pw_obs/np2;   t = ((1:np)*dt).*shp_obs(:,3)';

wr_obs= 2*pi*power_obs*shp_obs(:,2)/max(shp_obs(:,2));
wr_dec= 2*pi*power_dec*shp_dec(:,2)/max(shp_dec(:,2));

dphidt_obs = 2*pi*ramp_obs;
dphidt_dec = 2*pi*ramp_dec;

ph_obs = shp_obs(:,1)*pi/180 +  phase_obs + dphidt_obs*t' ;
ph_dec = shp_dec(:,1)*pi/180 +  phase_dec + dphidt_dec*t' ;

% Atualizar futuramente para np e tempo de duracao diferentes para o obs e
% o dec

U = 1;

for k=1:np 

    Hrf =   wr_obs(k)*(cos(ph_obs(k))*L_Mol(:,:,1) + sin(ph_obs(k))*L_Mol(:,:,2))...
          + wr_dec(k)*(cos(ph_dec(k))*L_Mol(:,:,3) + sin(ph_dec(k))*L_Mol(:,:,4));

    U = expm(-1i*(Hrf +  L_Mol(:,:,5) + L_Mol(:,:,6))*dt)*U;

end
[Ua, roh]= delay_Par(L_Spec(1)*1e-6,roha,L_Mol); 
roh = U*roh*U'; 
[Ub, roh]= delay_Par(L_Spec(2)*1e-6,roh,L_Mol);
end


