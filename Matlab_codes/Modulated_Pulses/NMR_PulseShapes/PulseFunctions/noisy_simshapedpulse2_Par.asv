function [U, roh] = noisy_simshapedpulse2_Par(shape_obs,shape_dec,pw_obs,pw_dec,... 
                    ramp_obs,ramp_dec,power_obs,power_dec,roha, L_Spec, L_Mol, noise)
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
% Both shapes (obs and dec) must have the same points numbers
% noise is a [p,q] = (1:5) vector where the parameters (px,py,pz, p,q), where pi's are pauli channel coefficients
%       and p,q are amplitude damping.

shp_obs = shape_obs;
shp_dec = shape_dec;

if isa(shape_obs,'char')
    shp_obs = load(shape_obs);
end
if isa(shape_dec,'char')
    shp_dec = load(shape_dec);
end

np = length(shp_obs(:,3));  np2 = sum(shp_obs(:,3));

dt = 1e-6*pw_obs/np2;   t = ((1:np)*dt).*shp_obs(:,3)';

wr_obs= 2*pi*power_obs*shp_obs(:,2)/max(shp_obs(:,2));
wr_dec= 2*pi*power_dec*shp_dec(:,2)/max(shp_dec(:,2));


dphidt_obs = 2*pi*ramp_obs;
dphidt_dec = 2*pi*ramp_dec;

ph_obs = shp_obs(:,1)*pi/180 + dphidt_obs*t' ;
ph_dec = shp_dec(:,1)*pi/180 + dphidt_dec*t' ;

% Atualizar futuramente para np e tempo de duracao diferentes para o obs e
% o dec

U = 1;

for k=1:np 

    Hrf =   wr_obs(k)*(cos(ph_obs(k))*L_Mol(:,:,3) + sin(ph_obs(k))*L_Mol(:,:,4))... 
          + wr_dec(k)*(cos(ph_dec(k))*L_Mol(:,:,5) + sin(ph_dec(k))*L_Mol(:,:,6));
            %wr(k)*(cos(ph(k))*mol.Ich{1,1} + sin(ph(k))*mol.Ich{1,2});
            %wr(k)*(cos(ph(k))*mol.Ich{2,1} + sin(ph(k))*mol.Ich{2,2});

    U = expm(-1i*(Hrf +  L_Mol(:,:,1) + L_Mol(:,:,2))*dt)*U;
    %U = expm(-1i*(Hrf +  mol.Hint + mol.Hzee)*dt)*U;

end
[Ua, roh]= delay_Par(L_Spec(1)*1e-6,roha,L_Mol); 
roh = U*roh*U'; 
[Ub, roh]= delay_Par(L_Spec(2)*1e-6,roh,L_Mol);

% apply noise
if ~isempty(noise)
    p = noise(1:3), q = noise(4:5);
    roh = Apply_composite_map(roh, p, q); %This won't work. I need to use then rho as the fidelity measure.
end

end


%-----------------%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%--------------------%

% Function to apply a combination of the three Pauli Maps + Amp. Damp.
% given the density matrix and the paramters of each map
% p(1) --> sx | p(2) --> sy | p(3) --> sz | p(4) and q --> Amp. Damp.
function rho_out = Apply_composite_map(rho_in, p, q)
    rho = rho_in;
    for j = 1:4
        KrausOps = KrausOps_gen(j, p(j), q);  % generate the Kaus operators for j-th map
        rho = Apply_Channel(rho, KrausOps);   % apply the j-th map
    end
    rho_out = rho;
end

% Function to apply a quantum channel given the density matrix and Kraus operators
function output = Apply_Channel(rho, krausOperators)
    output = zeros(size(rho));
    for i = 1:length(krausOperators)
        output = output + krausOperators{i} * rho * krausOperators{i}';
    end
end

% Generation of Kraus operators
function KrausOps = KrausOps_gen(channel, p, q)
    % Channel = 1 --> bit flip | 2 --> bit phase flip | 3 --> phase flip
    % Channel = 4 --> generalized amplitude damping | q = 0 | rho_inf -> |1>
    % q = 1 | rho_inf -> |0>
    % p = 1 full decoherence
    % q is related to temperature in the generalized amp. damp. when q in [0,1/2]
    id = [1.0 0.0; 0.0 1.0];
    s_k = {[0 1; 1 0] [0 -1i; 1i 0] [1 0; 0 -1] }; % {sx sy sz}
    % Kraus operators
    if channel == 1 || channel == 2 || channel == 3 % Pauli channels.
        M{1} = sqrt(1 - p/2)*id;
        M{2} = sqrt(p/2)*s_k{channel};
    end
    if channel == 4  % Gen amp damp
        M{1} = sqrt(q).*[1 0; 0 sqrt(1 - p)];
        M{2} = sqrt(q).*[0 sqrt(p); 0 0];
        M{3} = sqrt(1-q).*[sqrt(1 - p) 0; 0 1];
        M{4} = sqrt(1-q).*[0 0; sqrt(p) 0];
    end
    KrausOps = M;
end

function rho_out = ApplyOperator(rho_in, operator)
    rho_out = operator * rho_in * operator';
end

