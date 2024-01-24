clc;
clear all;
close all;

% Definicoes uteis
id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
ZZ = kron(sz,sz);

% projetores z
prj{1} = [1.0 0.0; 0.0 0.0];
prj{2} = [0.0 0.0; 0.0 1.0];

pw1 = 11.8;   % Tempo pulso de 90 observer (microsegundo)
pw2 = 10.5;   % Tempo pulso de 90 decopler (microsegundo)
JJ = 194.65;

addpath(['/home/user/nmr/Matlab_codes/nmr_lib/'])    % scropts para simulacao da evolucao e pulsos NMR
addpath(['/home/user/nmr/Matlab_codes/molecules/'])  % dados das moleculas
addpath(['/home/user/nmr/Matlab_codes/Math_util/'])  % dados das moleculas
%addpath([pwd '/home/user/nmr/Matlab_codes/Math_util/'])  % dados das moleculas
actdir =  pwd;

spectropar;
formato;
rho_a = prj{1};
p_0 = 0.7;
rho_b = [p_0 0.0; 0.0 (1-p_0)];
rho_ini = kron(rho_a,rho_b);

% Cnot otimizada - option 2
% + Ratacoes em Z na preaparaco do estado
[U1 rho] = rgpulse(pw1,90,0,rho_ini); % Rx(pi/2) x 1
[U2 rho] = delay((1/(2*JJ)),rho);  % Acoplamento J t = 3/(2J)
[U3 rho] = rgpulse(pw1,90,270,rho); % R-y(pi/2) x 1
%U_part1 = U3*U2*U1
