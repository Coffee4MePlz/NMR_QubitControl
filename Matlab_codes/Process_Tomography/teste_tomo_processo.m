clc
clear variables
close all

% Start - adding path to useful scripts
add_scripts(pwd) %here pwd should be directory ./nmr/

id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
ZZ = kron(sz,sz);

syms thet1 thet2
Rotx(thet1) = [cos(thet1/2) -1i*sin(thet1/2); -1i*sin(thet1/2) cos(thet1/2)];
Roty(thet2) = [cos(thet2/2) -sin(thet2/2); sin(thet2/2) cos(thet2/2)];

op = {1i*id sx sy sz};
% Defining the Input states MUB 
input_states = {[1 0; 0 0] [0 0; 0 1] [1/2 1/2; 1/2 1/2] [1/2 -1i/2; 1i/2 1/2]}; % |0> |1> |x+> |y->
% Output experimental states

%% Unitary evolution
% operator = Rotx(pi/2);
% output_states = {ApplyOperator(input_states{1}, operator) ApplyOperator(input_states{2}, operator)...
%   ApplyOperator(input_states{3}, operator) ApplyOperator(input_states{4}, operator)};



% Non-unitary evolution
% Decoherence_Channels(rho_in, p, q, channel)
% Channel = 1 --> bit flip | 2 --> bit phase flip | 3 --> phase flip | 4 --> depolarizing
% Channel = 5 --> amplitude damping | 6 --> gneralized amplitude damping | 7 --> phase damping 
% p = 0.7; q = 0.0; channel = 5;
% output_states = {Decoherence_Channels(input_states{1}, p, q, channel) ...
%     Decoherence_Channels(input_states{2}, p, q, channel) Decoherence_Channels(input_states{3}, p, q, channel) ...
%     Decoherence_Channels(input_states{4}, p, q, channel)};

% Unitary + Non-unitary evolution
% Decoherence_Channels(rho_in, p, q, channel)
% Channel = 1 --> bit flip | 2 --> bit phase flip | 3 --> phase flip | 4 --> depolarizing
% Channel = 5 --> amplitude damping | 6 --> generalized amplitude damping | 7 --> phase damping 
% 
operator = Rotx(pi/2); % parte unitaria
p = 0.03; q = 0.03; channel = 1;

output_states = {Decoherence_Channels(ApplyOperator(input_states{1}, operator), p, q, channel) ...
     Decoherence_Channels(ApplyOperator(input_states{2}, operator), p, q, channel) ...
     Decoherence_Channels(ApplyOperator(input_states{3}, operator), p, q, channel) ...
     Decoherence_Channels(ApplyOperator(input_states{4}, operator), p, q, channel)};
% load(['/home/quantum/nmr/Matlab_codes/Carlos/rho_tomo_proc'],'rho_1q_00','rho_1q_10','rho_1q_xplus',...
%     'rho_1q_yplus','rho_1q_Rot00','rho_1q_Rot10','rho_1q_Rotxplus','rho_1q_Rotyplus');

% input_states = {rho_1q_00 rho_1q_10{1} rho_1q_xplus{1} rho_1q_yplus{1}};
% output_states = {rho_1q_Rot00{1} rho_1q_Rot10{1} rho_1q_Rotxplus{1} rho_1q_Rotyplus{1}};
 

% Generic Choi Matrix
syms c11 c12 c13 c14 c21 c22 c23 c24 c31 c32 c33 c34 c41 c42 c43 c44
choi_matrix_sym = [c11 c12 c13 c14; c21 c22 c23 c24; c31 c32 c33 c34; c41 c42 c43 c44];

% % For didatics only
% for k = 1:4
%   result_1{k} = zeros(2);
%   for i = 1:4
%     for j = 1:4
%       result_1{k} = result_1{k} + choi_matrix_sym(i,j)*op{i}*input_states{k}*ctranspose(op{j});
%     end
%   end
% end


% Considering 1_qubit generic input
syms r11 r12 r21 r22
input_gen = [r11 r12; r21 r22];

% Generic process tranfomation
result_gen = zeros(2);
for i = 1:4
  for j = 1:4
    result_gen = result_gen + choi_matrix_sym(i,j)*op{i}*input_gen*ctranspose(op{j});
  end
end

%  Coeficient matrix Sum_i,j Xi_i,j G(i) * rho_in * G(J) == rho_out 
syms r_out11 r_out12 r_out21 r_out22
eqn1 = result_gen(1,1) == r_out11;
eqn2 = result_gen(1,2) == r_out12;
eqn3 = result_gen(2,1) == r_out21;
eqn4 = result_gen(2,2) == r_out22;
% Puting in the matrix form 4x16 for a given input state
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], ...
  [c11, c12, c13, c14, c21, c22, c23, c24, c31, c32, c33, c34, c41, c42, c43, c44]);
% The matrix form for the 4 MUB states
for k = 1:4
  A_aux{k} = double(subs(A, {r11, r12, r21, r22}, {input_states{k}(1,1),...
    input_states{k}(1,2), input_states{k}(2,1), input_states{k}(2,2)}));
  B_aux{k} = vpa(subs(B, {r_out11, r_out12, r_out21, r_out22},...
    {output_states{k}(1,1), output_states{k}(1,2), output_states{k}(2,1), output_states{k}(2,2)}));
end

%
% Concatenating the 4 evolutions for the 4 initial state resulting in the system A * C = B now as 16 x 16 matriiiix
% C is the linearized Choi matrix as C = [c11 c12 c13 c14 c21 c22 c23 c24 c31 c32 c33 c34 c41 c42 c43 c44]' 
A_conc = [A_aux{1}; A_aux{2}; A_aux{3}; A_aux{4}];
B_conc = [B_aux{1}; B_aux{2}; B_aux{3}; B_aux{4}];

Choi_coef = linsolve(A_conc,B_conc);

% Choi (matrix) for the map representation as rho_out = Sum_i,j chi_i,j G_i * rho_in * G_j'

Choi_matrix = round([Choi_coef(1) Choi_coef(2) Choi_coef(3) Choi_coef(4);...
  Choi_coef(5) Choi_coef(6) Choi_coef(7) Choi_coef(8);...
  Choi_coef(9) Choi_coef(10) Choi_coef(11) Choi_coef(12);...
  Choi_coef(13) Choi_coef(14) Choi_coef(15) Choi_coef(16)],4)

[vv lambda] = eig(Choi_matrix)

% Operator Sum representation (Kraus)  rho_out = Sum_i lambda_i * L_i * rho_in * G_i'
% lambda are the eigevalues of the Choi matrix chi
%% Kraus operators KK(j) = lambda_i * L_i
Soma_Kraus = zeros(2);
for j = 1:4
    if lambda(j,j) >= 0.000001
        ll(j) = lambda(j,j);
        Kraus{j} = (vv(1,j)*op{1} + vv(2,j)*op{2} + vv(3,j)*op{3} + vv(4,j)*op{4});
        Soma_Kraus = Soma_Kraus + lambda(j,j)*ctranspose(Kraus{j})*Kraus{j};  
    end
end
% Soma_Kraus should be = eye(2) 
round(Soma_Kraus,2)



%% Testes
PlotProcess(Choi_matrix,1)

% 
% rho_out_map = Apply_Map(input_states{3},Choi_matrix);
% 
% round(rho_out_map,4)

%Unitaria = vec(1,1)*op{1} + vec(2,1)*op{2} + vec(3,1)*op{3} + vec(4,1)*op{4}


% --- Functions

function result = add_scripts(str)
str1 = {};
    str1{1} = append(str,'/Matlab_codes/nmr_lib/');
    str1{2} = append(str,'/Matlab_codes/gates/');
    str1{3} = append(str, '/Matlab_codes/Math_util/');
    for i =1:length(str1)
        addpath([str1{i}])
    end
    %addpath(['/home/quantum/nmr/Matlab_codes/nmr_lib/'])    % scripts para simulacao da evolucao e pulsos NMR
    %addpath(['/home/quantum/nmr/Matlab_codes/gates/'])      % gates
    %addpath(['/home/quantum/nmr/Matlab_codes/Math_util/'])  % math utils 
end

function rho_out = ApplyOperator(rho_in, operator)
rho_out = operator * rho_in * operator';
end

function rho_out = Decoherence_Channels(rho_in, p, q, channel)
% Channel = 1 --> bit flip | 2 --> bit phase flip | 3 --> phase flip | 4 --> depolarizing
% Channel = 5 --> amplitude damping | 6 --> gneralized amplitude damping | 7 --> phase damping 
% p = 1 fill decoherence
% q is related to temp in generalized amp. damp.
id = [1.0 0.0; 0.0 1.0];
%sz = [1.0 0.0; 0.0 -1.0]; sx = [0 1; 1 0]; sy = [0 -1i; 1i 0];
s_k = {[0 1; 1 0] [0 -1i; 1i 0] [1 0; 0 -1] }; % {sx sy sz}

% Kraus operators
if channel == 1 || channel == 2 || channel == 3 % Pauli chan.
  M{1} = sqrt(1 - p/2).*id;
  M{2} = sqrt(p/2).*s_k{channel};
end
if channel == 4                 % depolarazing
  M{1} = sqrt(1 - 3*p/4).*id;
  M{2} = sqrt(p/4).*s_k{1};
  M{3} = sqrt(p/4).*s_k{2};
  M{4} = sqrt(p/4).*s_k{3};
end
if channel == 5
  M{1} = [1 0; 0 sqrt(1 - p)];
  M{2} = [0 sqrt(p); 0 0];result_gen = zeros(2);
end
if channel == 6
  M{1} = sqrt(q).*[1 0; 0 sqrt(1 - p)];
  M{2} = sqrt(q).*[0 sqrt(p); 0 0];
  M{3} = sqrt(1-q).*[sqrt(1 - p) 0; 0 1];
  M{4} = sqrt(1-q).*[0 0; sqrt(p) 0];
end
if channel == 7
  M{1} = [1 0; 0 sqrt(1 - p)];
  M{2} = [0 0; 0 sqrt(p)];
end

rho_out = zeros(2);
for k = 1:length(M)
  rho_out = rho_out + M{k}*rho_in*ctranspose(M{k});
end

end

function rho_out = Apply_Map(rho_in,Choi_matrix)
id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];

op = {id sx sy sz};

result = zeros(2);
for i = 1:4
  for j = 1:4
    result = result + Choi_matrix(i,j)*op{i}*rho_in*ctranspose(op{j});
  end
end

rho_out = result

end
