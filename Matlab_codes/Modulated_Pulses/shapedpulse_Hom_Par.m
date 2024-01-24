function [U, roh] = shapedpulse_Hom_Par(shape_obs,pw_obs,phase_obs,ramp_obs,power_obs,roha, L_Spec, L_Mol)
%(shape_obs,shape_dec,pw_obs,pw_dec,phase_obs,phase_dec,... 
    %ramp_obs,ramp_dec,power_obs,power_dec,roha) % future inputs for when developing for dec and obs
% for now this is only on obs, but can be changed to have dec as well, like as in
% [U, roh] = simshapedpulse(shape_obs,shape_dec,pw_obs,pw_dec,phase_obs,phase_dec,ramp_obs,ramp_dec,power_obs,power_dec,roha)
% U ->  Pulse propagator (this incluies the free evoluiton induced by the rofs times)
% Both shapes (obs and dec) must have the same points numbers (for now)
% L_... variables used in parallelization
arguments
    shape_obs (:,3) 
    pw_obs double = 200 %pulse width in microseconds
    phase_obs (:,1) = [] %phase of the pulse in degree, default is []
    ramp_obs double = 0 %introduce a linear ramp on the phase to shift the excitation region (given in Hz)
    power_obs double = 20*1e3 %it is the maximum amplitude of the pulse (B1 in Hz) 
    roha (:,:) = [1 0; 0 0] %density matrix before the pulse, after the pulse its `roh`
    L_Spec (1,:) =  [10.0000   24.8750] 
    L_Mol (:,:,:) = zeros(8,8,6)
end

shp_obs = shape_obs;

if isa(shape_obs,'char')
    shp_obs = load(shape_obs);
end

np = length(shp_obs(:,3));  np2 = sum(shp_obs(:,3));

dt = 1e-6*pw_obs/np2;   t = ((1:np)*dt).*shp_obs(:,3)';

wr_obs= 2*pi*power_obs*shp_obs(:,2)/max(shp_obs(:,2)); 
%wr_dec= 2*pi*power_dec*shp_dec(:,2)/max(shp_dec(:,2)); %not used in single pulse

dphidt_obs = 2*pi*ramp_obs;
%dphidt_dec = 2*pi*ramp_dec; %not used in single pulse

% code for overriding phase from shape matrix when phase is explicity inputed
if (~isempty(phase_obs))
    ph_obs = phase_obs*pi/180 +dphidt_obs*t';
else
    ph_obs = shp_obs(:,1)*pi/180 + dphidt_obs*t' ; 
end

%ph_dec = shp_dec(:,1)*pi/180 +  phase_dec + dphidt_dec*t' ; %not used in single pulse

U = 1;
for k=1:np 
Hrf= 0;

   [p1,p2,p3] = size(L_Mol);
   Ns = (p3-2)/2;
   for j=1:Ns
        Hrf = Hrf + wr_obs(k)*(cos(ph_obs(k))*L_Mol(:,:,(2*j+1)) + sin(ph_obs(k))*L_Mol(:,:,(2*j+2)));
   end
   
   U = expm(-1i*(Hrf +  L_Mol(:,:,1) + L_Mol(:,:,2))*dt)*U; % mol.Hint + mol.Hzee)*dt)*U;

end
%When unparallel
%{
[Ua roh]= delay(spectro.rof1*1e-6,roha); 
roh = U*roh*U'; 
[Ub roh]= delay(spectro.rof2*1e-6,roh);
%}

%when Parallel :
[Ua, roh]= delay_Par(L_Spec(1)*1e-6,roha,L_Mol); 
roh = U*roh*U'; 
[Ub, roh]= delay_Par(L_Spec(2)*1e-6,roh,L_Mol);

end


