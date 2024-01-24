clc
clear variables
close all

% Start - adiciona os paths do scripts importantes
start;
% Definicoes uteis
%%
id = eye(2);
sp = {[0 1; 1 0] [0 -1i; 1i 0] [1.0 0.0; 0.0 -1.0]}; % Pauli matrix
Rot{1} = @(thet) [cos(thet/2) -1i*sin(thet/2); -1i*sin(thet/2) cos(thet/2)]; % Rx
Rot{2} = @(thet) [cos(thet/2) -sin(thet/2); sin(thet/2) cos(thet/2)]; % Ry
Rot{3} = @(thet) [exp(-1i*thet) 0; 0 exp(1i*thet)]; % Rz
%% simulation test
% rho_in = 1/2*[1 1i; -1i 1]; %estado inicial 
% p = [1/2, 0, 0, 0]; % Parametros px py pz (Pauli maps) + p amp. damp.
% q = 1; % Parametro q amp. damp. 
% rho_out = Apply_composite_map(rho_in, p, q) % combined maps

%defining parameters:
% true par = [global phase, \theta_x, \theta_y, \theta_z, p_x, p_y, p_z, p_a, q_a ] 
true_par = [0 0 pi/2 0 0.0 0.1 0 0.0 0.0]

% Tomography input
% Defining the Input states MUB 
input_states = {[1 0; 0 0] [0 0; 0 1] [1/2 1/2; 1/2 1/2] [1/2 1i/2; -1i/2 1/2]}; % |0> |1> |x+> |y->

% Output teo/experimental states
% Unitary evolution
thet_global = true_par(1); thet1 = true_par(2);thet2 =true_par(3); thet3 = true_par(4);


% plotting to check
%rhotoplot{1} = output_states{4};
%PlotBloch(rhotoplot,1,1, [1,0,1], false)
%

% general unitary: 
unitary = exp(-1i*thet_global)*Rot{3}(thet3)*Rot{2}(thet2)*Rot{1}(thet1); 
% example: hadamard
%unitary =  1/sqrt(2)*[1 1 ; 1 -1];

for i = 1:length(input_states)
    output_states{i} = ApplyOperator(input_states{i}, unitary); 
end


% Non-unitary evolution
p = [true_par(5), true_par(6), true_par(7), true_par(8)]; % Parametros px py pz (Pauli maps) + p(4) amp. damp.
q = true_par(9); % Parametro q amp. damp. 
for i = 1:length(input_states)
    output_states{i} = Apply_composite_map(output_states{i}, p, q);
end

[Choi_matrix,Kraus,ll,Soma_Kraus] = process_tomography(input_states,output_states);
%%
%PlotProcess2(Choi_matrix,1)

Choi_matrix_exp = Choi_matrix;
%par = ones(1,9);
par = zeros(1,9);
%par(2) = pi;
lb = [-pi -pi -pi -pi 0 0 0 0 0];%- 0.001;
ub = [pi pi pi pi 0.1 0.1 0.1 0.1 0.1] + 0.001;

f = @(par) chi_proc(Choi_matrix_exp, par);
m=1; Mit =100; %chi = ones(1,Mit); par_opt = cell(1,Mit); par_ini = cell(1,Mit);
for m =1:Mit
    par_opt{m} = zeros(1,length(par));
    par_ini{m} = zeros(1,length(par));
end
tic
parfor m = 1:Mit
    %par_ini{m} = zeros(1,length(par));
    par_ini{m} = set_par_ini(Choi_matrix_exp, lb,ub,m);
 
    % Set optimization options
    options = optimset('MaxFunEvals', 1e6,'MaxIter', 1e6,...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-7);
    [par_opt{m},chi(m),~,~] = fminsearchbnd(f,par_ini{m},lb,ub,options);
end

toc
[mm ind] = min(chi);
par_opt{ind}
%%
Choi_mat = find_process(par_opt{ind});
PlotProcess2(Choi_mat,3)
PlotProcess2(Choi_matrix_exp,1)
%% EXPERIMENTAL TESTTS

% Tomography input
% Defining the Input states MUB 
input_states = {[1 0; 0 0] [0 0; 0 1] [1/2 1/2; 1/2 1/2] [1/2 -1i/2; 1i/2 1/2]}; % |0> |1> |x+> |y->

% Output teo/experimental states

% EXPERIMENTAL DATA:
%%{
% 1) fectch experimental data at teste_tomo_processo.m
%exp_output = rho_1q_Rotymod2;
%exp_output = rho_ml_full;
exp_output = rhoTH_fullmod;
% 2) writing input and output states explicitly:
%input_states = input_states;
%output_states = cell(1,length(exp_output))
%for k =  1:11
for j = 1:4
    input_states{j} = exp_output{j,1};
    output_states{j} = exp_output{j,k};
end
%}

[Choi_matrix,Kraus,ll,Soma_Kraus] = process_tomography(input_states,output_states);
PlotProcess2(Choi_matrix,1)
Choi_matrix_exp = Choi_matrix;


par = zeros(1,9);
lb = [-pi -pi -pi -pi 0 0 0 0 0];%- 0.001;
ub = [pi pi pi pi 0.1 0.1 0.1 0.1 0.1] + 0.001;

f = @(par) chi_proc(Choi_matrix_exp, par);
m=1; Mit =200; %chi = ones(1,Mit); par_opt = cell(1,Mit); par_ini = cell(1,Mit);
parfor m =1:Mit
    par_opt{m} = zeros(1,length(par));
    par_ini{m} = zeros(1,length(par));
end
tic
parfor m = 1:Mit
    par_ini{m} = set_par_ini(Choi_matrix_exp, lb,ub,m);
 
    % Set optimization options
    options = optimset('MaxFunEvals', 1e6,'MaxIter', 1e6,...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-7);
    [par_opt{m},chi(m),~,~] = fminsearchbnd(f,par_ini{m},lb,ub,options);
end

toc
[mm ind] = min(chi);
par_opt{ind}

%%
Choi_mat = find_process(par_opt{ind});
PlotProcess2(Choi_mat,3)
PlotProcess2(Choi_matrix_exp,1)
%M_ChoiMat{k} = Choi_matrix_exp;
%optimal_par{k} = par_opt{ind};
%disp(k)
%end

%% some plotting routines

for j=1:11
    PauliVectors(j,:) = optimal_par{j}(5:9);
end
t = 0:20:201;

figure(5); % Opens a new figure window
hold on; % Allows multiple plots in the same figure
for i = 1:size(PauliVectors, 2)
    if i <= 3
        plot(t, PauliVectors(:, i), '-o', 'DisplayName', sprintf('p_%d', i), 'LineWidth', 3);
    elseif i == 4
        plot(t, PauliVectors(:, i), '-o', 'DisplayName', 'p_a','LineWidth', 3 );
    else
        plot(t, PauliVectors(:, i), '-o', 'DisplayName', 'q_a', 'LineWidth', 3);
    end
end

hold off;
grid on
% Add labels and legend
xlabel('Time (us)');
ylabel('Value');
ylim([-0.01,0.15]);
title('Pauli Noise Vectors and amplitude damping Over Time');
legend('show'); % Displays a legend
%%
PlotProcess2(M_ChoiMat{11},8)



%% GENERATING DATASET
%defining parameters:
lb = [-pi -pi -pi -pi 0 0 0 0 0];%- 0.001;
ub = [pi pi pi pi 0.2 0.2 0.2 0.2 0.2] + 0.001;
path_matrix = "/home/user/nmr/Matlab_codes/Process_Tomography/ML-QPT-Dataset/QPT_Matrices.csv";
path_labels = "/home/user/nmr/Matlab_codes/Process_Tomography/ML-QPT-Dataset/QPT_labels.csv";
for ii = 2:200001
    true_par = random_true_par(lb,ub);
    %unitary case:
    %true_par(5:9) = zeros(1,5);
    
    % Tomography input
    % Defining the Input states MUB 
    input_states = {[1 0; 0 0] [0 0; 0 1] [1/2 1/2; 1/2 1/2] [1/2 -1i/2; 1i/2 1/2]}; % |0> |1> |x+> |y->
    
    % Output teo/experimental states
    % Unitary evolution
    thet_global = true_par(1); thet1 = true_par(2);thet2 =true_par(3); thet3 = true_par(4); 
    unitary = exp(-1i*thet_global)*Rot{3}(thet3)*Rot{2}(thet2)*Rot{1}(thet1); 
    for i = 1:length(input_states)
        output_states{i} = ApplyOperator(input_states{i}, unitary); 
    end
    
    % Non-unitary evolution
    p = [true_par(5), true_par(6), true_par(7), true_par(8)]; % Parametros px py pz (Pauli maps) + p(4) amp. damp.
    q = true_par(9); % Parametro q amp. damp. 
    for i = 1:length(input_states)
        output_states{i} = Apply_composite_map(output_states{i}, p, q);
    end
    [Choi_matrix,~,~,~] = process_tomography(input_states,output_states);
    

    if ii==1
        dlmwrite(path_matrix, Choi_matrix);
        dlmwrite(path_labels, true_par);
    else
        dlmwrite(path_matrix, Choi_matrix, '-append');
        dlmwrite(path_labels, true_par, '-append');
    end
    if mod(ii,10000) == 0
        disp(ii)
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>>>>>>>>>>>>>>>>>>>>> functions <<<<<<<<<<<<<<<<<<<<<<<<<<

function true_par = random_true_par(lb,ub)
    for k =1:length(lb)
        rn = rand(1,9);
        true_par = lb+ (ub-lb).*rn;
    end
end

function par = set_par_ini(Choi_matrix_exp, lb, ub, m)
    par = zeros(1,9);
    preff_dir = find(diag(Choi_matrix_exp) > 0);
    if m ==1
        for i=1:length(preff_dir)
            angle = Choi_matrix_exp(preff_dir(i),preff_dir(i))*pi;
            if  angle>0
                par(preff_dir(i)) = angle;
            end
        end
    else
        for i=1:length(preff_dir)
            angle = Choi_matrix_exp(preff_dir(i),preff_dir(i))*pi;
            if  angle>0
                par(preff_dir(i)) = angle*(1-normrnd(0,1/4));
            end
        end
        for j=5:9
            rnd1 = rand(1);
            par(j) = rnd1*lb(j) + (1-rnd1)*ub(j);
        end
    end
end

function chi = chi_proc(Choi_matrix_exp, par)

    Choi_mat = find_process(par);
    %anglesum = 0;
    %langle = find(diag(Choi_matrix_exp) > 0);
    
    % think about minimizing with greater strenght the angle associated to biggest element in choi matrix
    chi = sum(sum( real(Choi_mat-Choi_matrix_exp).^2 + imag(Choi_mat-Choi_matrix_exp).^2 ));% + anglesum;

end

function Choi_mat = find_process(par)

    id = eye(2);
    sp = {[0 1; 1 0] [0 -1i; 1i 0] [1.0 0.0; 0.0 -1.0]}; % Pauli matrix
    Rot{1} = @(thet) [cos(thet/2) -1i*sin(thet/2); -1i*sin(thet/2) cos(thet/2)]; % Rx
    Rot{2} = @(thet) [cos(thet/2) -sin(thet/2); sin(thet/2) cos(thet/2)]; % Ry
    Rot{3} = @(thet) [exp(-1i*thet) 0; 0 exp(1i*thet)]; % Rz
    % Tomography input
    % Defining the Input states MUB 
    input_states = {[1 0; 0 0] [0 0; 0 1] [1/2 1/2; 1/2 1/2] [1/2 -1i/2; 1i/2 1/2]}; % |0> |1> |x+> |y->
    % Output teo/experimental states
    % Unitary evolution
    thet_global = par(1); thet1 = par(2);thet2 =par(3); thet3 = par(4);
    unitary = exp(-1i*thet_global)*Rot{3}(thet3)*Rot{2}(thet2)*Rot{1}(thet1); 
    for i = 1:length(input_states)
        output_states{i} = ApplyOperator(input_states{i}, unitary); 
    end
    % Non-unitary evolution
    p = [par(5), par(6), par(7), par(8)]; % Parametros px py pz (Pauli maps) + p(4) amp. damp.
    q = par(9); % Parametro q amp. damp. 
    for i = 1:length(input_states)
        output_states{i} = Apply_composite_map(output_states{i}, p, q);
    end
    
    [Choi_mat,Kraus,ll,Soma_Kraus] = process_tomography(input_states,output_states);

end



% Process Tomography
function [Choi_matrix,Kraus,ll,Soma_Kraus] = process_tomography(input_states,output_states)
    % Process tomography for the operator basis 
    op = {1i*eye(2) [0.0 1.0; 1.0 0.0] [0.0 -1i; 1i 0.0] [1.0 0.0; 0.0 -1.0]};
    
    % A*Choi_coef = B
    syms r11 r12 r21 r22 r_out11 r_out12 r_out21 r_out22
    % A = [r11, r12,  r12*1i,  r11, r21, r22,  r22*1i,  r21, -r21*1i, -r22*1i,  r22, -r21*1i,  r11,  r12,  r12*1i,  r11;...
    %     r12, r11, -r11*1i, -r12, r22, r21, -r21*1i, -r22, -r22*1i, -r21*1i, -r21,  r22*1i,  r12,  r11, -r11*1i, -r12;...
    %     r21, r22,  r22*1i,  r21, r11, r12,  r12*1i,  r11,  r11*1i,  r12*1i, -r12,  r11*1i, -r21, -r22, -r22*1i, -r21;...
    %     r22, r21, -r21*1i, -r22, r12, r11, -r11*1i, -r12,  r12*1i,  r11*1i,  r11, -r12*1i, -r22, -r21,  r21*1i,  r22];
    A = [r11, r12*1i, -r12,  r11*1i, -r21*1i, r22,  r22*1i,  r21, -r21, -r22*1i,  r22, -r21*1i, -r11*1i,  r12,  r12*1i,  r11;...
        r12, r11*1i,  r11, -r12*1i, -r22*1i, r21, -r21*1i, -r22, -r22, -r21*1i, -r21,  r22*1i, -r12*1i,  r11, -r11*1i, -r12;...
        r21, r22*1i, -r22,  r21*1i, -r11*1i, r12,  r12*1i,  r11,  r11,  r12*1i, -r12,  r11*1i,  r21*1i, -r22, -r22*1i, -r21;...
        r22, r21*1i,  r21, -r22*1i, -r12*1i, r11, -r11*1i, -r12,  r12,  r11*1i,  r11, -r12*1i,  r22*1i, -r21,  r21*1i,  r22];
     
    
    B = [r_out11; r_out12; r_out21; r_out22];
    
    for k = 1:4
      A_aux{k} = double(subs(A, {r11, r12, r21, r22}, {input_states{k}(1,1),...
        input_states{k}(1,2), input_states{k}(2,1), input_states{k}(2,2)}));
      B_aux{k} = double(subs(B, {r_out11, r_out12, r_out21, r_out22},...
        {output_states{k}(1,1), output_states{k}(1,2), output_states{k}(2,1), output_states{k}(2,2)}));
    end
    
    % Concatenating the 4 evolutions for the 4 initial state resulting in the system A * C = B now as 16 x 16 matriiiix
    % Now C is the linearized Choi matrix as C = [c11 c12 c13 c14 c21 c22 c23 c24 c31 c32 c33 c34 c41 c42 c43 c44]' 
    A_conc = [A_aux{1}; A_aux{2}; A_aux{3}; A_aux{4}];
    B_conc = [B_aux{1}; B_aux{2}; B_aux{3}; B_aux{4}];
    
    Choi_coef = linsolve(A_conc,B_conc);
    
    % Choi (matrix) for the map representation as rho_out = Sum_i,j chi_i,j G_i * rho_in * G_j'
    Choi_matrix = round([Choi_coef(1) Choi_coef(2) Choi_coef(3) Choi_coef(4);...
      Choi_coef(5) Choi_coef(6) Choi_coef(7) Choi_coef(8);...
      Choi_coef(9) Choi_coef(10) Choi_coef(11) Choi_coef(12);...
      Choi_coef(13) Choi_coef(14) Choi_coef(15) Choi_coef(16)],4);
    
    [vv lambda] = eig(Choi_matrix); % Eigenvec and eigenvals of the Choi matrix
    
    % Operator Sum representation (Kraus)  rho_out = Sum_i lambda_i * L_i * rho_in * G_i'
    % lambda(j,j) are the eigevalues of the Choi matrix chi
    %% Kraus operators KK(j) = lambda_i * L_i
    Soma_Kraus = zeros(2);
    for j = 1:4
        if lambda(j,j) >= 0.0000001
            ll(j) = lambda(j,j);
            Kraus{j} = sqrt(ll(j))*(vv(1,j)*op{1} + vv(2,j)*op{2} + vv(3,j)*op{3} + vv(4,j)*op{4});
            Soma_Kraus = Soma_Kraus + ctranspose(Kraus{j})*Kraus{j};  
        end
    end
    % Soma_Kraus should be = eye(2) 
end


% -----------------------------
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

% -----------------------------
% Function to apply a quantum channel given the density matrix and Kraus operators
function output = Apply_Channel(rho, krausOperators)
    output = zeros(size(rho));
    for i = 1:length(krausOperators)
        output = output + krausOperators{i} * rho * krausOperators{i}';
    end
end

% ----------------------------
% Genration of Kraus operators
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

% ----------------------------------
function rho_out = ApplyOperator(rho_in, operator)
    rho_out = operator * rho_in * operator';
end

