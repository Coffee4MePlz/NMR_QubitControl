%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TUTORIAL ON HOW TO USE MODULATED PULSES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATE 23/01/2024 %%%
%{ % -->> %
% In this section we will give some worked examples on how to use the file on modulated pulses
%           modulados_HomQbits_Par_withNoise.m
% one can just copy and paste them before calling the functions 
%
% Tutorial Hints are labeled with "% -->> %"
%}


% % -->> % as nice practice, one can clear the memory with:
% clearing memory
clc;
clear all;
close all;

% % -->> % Choose a molecule (can be cloroformio, trifluor1,2,3, formato, etc) and set the spectropar variable
% fetching objects
spectropar;
cloroformio;
%trifluor3;
%trifluor2;
%trifluor1;
%formato
%tricloro;

% % -->> % sets the global variables mol and spectro to local ones, L_Spec,L_Mol
global mol spectro
[L_Spec, L_Mol] = fetch_global_val(spectro, mol);

% % -->> % Declaring some useful variables
% Useful matrices % 
id = [1 0; 0 1];
sz = [1 0; 0 -1];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
ZZ = kron(sz,sz);

% projetores z
prj{1} = [1.0 0.0; 0.0 0.0];
prj{2} = [0.0 0.0; 0.0 1.0];
% projetores x

prjx{1} = 0.5*[1.0 1.0; 1.0 1.0];      %|+><+|
prjx{2} = 0.5*[1.0 -1.0; -1.0 1.0];    %|-><-|
% projetores y
prjy{1} = 0.5*[1.0 -1i; 1i 1.0];      %|+><+|
prjy{2} = 0.5*[1.0 1i; -1i 1.0];    %|-><-|



%% EXAMPLE 1 : modulating for a C-NOT decomposed unitaries in Heteronuclear systems


%fetching molecule properties
cloroformio
global mol spectro
[L_Spec, L_Mol] = fetch_global_val(spectro,mol);

% % -->> % NOW we have to set the pulse parameters, It will go to a huge structure,
% % -->> % but can be set first here if there is much redundancy

%defining pulse variables
pulse_width = 200 ; np=200;
%experimental maxpowers: Can be set to half
maxpower_obs = 21030; %1/2*21030; %according to hamiltonian estimation 
maxpower_dec = 23474; %1/2*23474; % 
minpower = 0*1e3;
power_bound_obs = [minpower maxpower_obs]; 
power_bound_dec = [minpower maxpower_dec];
acc = 1e-4; %numerical accuracy (its an optimization parameter)
s_a = 4; % Number of coefficients for amplitude Fourier Series
s_p = 14; % Number of coefficients for amplitude Fourier Series
Par = [2 200]; % Defining the type of parallelization, "[0 10]" means no paralelization, "[1 10]" means 
% parallelization defined outside the function as a parfor (paralelized for loop) with 10 paralel runs, and "[2 10]" parallelization
% inside as a Multistart() MATLAB function, with 10 paralel runs.,
verbose = 1; %verbose is 0 if you don't want to print anything, and 1 to print development

% % -->> % WE can set an initial value for the density matrix
rhoini = kron(prj{1},prj{1});
rho_evol{1} = rhoini;

% % -->> % WE will now define the Circuit we wish to run

% CNOT FOR NMR has the following circuit
% % -->> % DEFINE the target unitaries

%U0 = kron(id,expm(-1i*pi/2*sx)); % OPTIONAL
U1 = kron(id, expm(-1i*pi/4*sy)); % flip Second qubit
U2 = expm(-1i*pi/4*kron(sz,sz)); % interaction gate
U3 = kron(id,expm(1i*pi/4*sx)); % all other gates after interaction
U4 = kron(expm(1i*pi/4*sy),expm(1i*pi/4*sy));
U5 = kron(expm(1i*pi/4*sx),expm(-1i*pi/4*sx));
U6 = kron(expm(-1i*pi/4*sy),expm(-1i*pi/4*sy));

U = U6*U5*U4*U3*U2*U1; % global unitary
% % -->> % The circuit is set to be done in three gates, we this is given by the ordering of QbitOrder, which is a cell
% % -->> % {"sim"} is for simultaneous pulses (both channels), {"int"} is interaction gate
% % -->> % {"sing",i} is for single qubit pulsing, and i=0,1,2,3,.. is the target qubit
% % -->> %  {"grad"} will apply a gradient, i.e., delete the non-diagonal terms of \rho matrix
QbitOrder = {{"sim"},{"int"},{"sim"}}; 
% storting unitaries as as (U, n) 3D matrix (an ordered tensor), following the circuit order in modulated pulses
U_target(:,:,1) = U1;
U_target(:,:,2) = U2;
U_target(:,:,3) = U6*U5*U4*U3;


%%%%%%%%%%%%  For state preparation the lines above can be substituted by the following batch of code %%%%%%%%%
% % -->> % ANOTHER example with state preparation
% state preparation: following Batalhão page 107.
%{
% free gates
U_free2 = expm(-1i*(pi/4)*kron(sz,sz)); %U_free2 = expm(-1i*(3*pi/4)*kron(sz,sz));
U_free1 = expm(-1i*(pi/8)*kron(sz,sz)); %U_free1 = expm(-1i*(3*pi/8)*kron(sz,sz));
rhoini = 4*kron(sz,id) + kron(id,sz);
rho_evol{1} =rhoini;

m =0; n =0;
k = (-1)^(m); j = (-1)^(n);
% defining Unitary order:
QbitOrder = {{"sim",1},{"int"},{"sim",1},{"int"},{"sim",1},{'grad'},{"sim",1},{"int"},{"sim",1},{'grad'},{"sim",1}};
    U_target(:,:,1) = kron(expm(1i*pi/4*sx), expm(1i*pi/4*sx));
    U_target(:,:,2) = U_free1 ;
    U_target(:,:,3) = kron(expm(-1i*pi/4*sy), expm(-1i*pi/4*sy));
    U_target(:,:,4) = U_free1;
    U_target(:,:,5) = kron(expm(-k*1i*pi/4*sx), expm(-j*1i*pi/4*sx));
    U_target(:,:,6) = eye(4); %gradient
    U_target(:,:,7) = kron(expm(1i*pi/8*sy),expm(1i*pi/8*sy));
    U_target(:,:,8) = U_free1;
    U_target(:,:,9) = kron(expm(j*1i*pi/12*sx),expm(k*1i*pi/12*sx));
    U_target(:,:,10) = eye(4); %gradient
    U_target(:,:,11) = kron(expm(-1i*pi/8*sx),id);
%vizualization of evolution
U =1;
for j =1:length(QbitOrder)
    currentQbit = QbitOrder{j};
    if (currentQbit == 'grad')
        D = diag(rho_evol{j});
        rho1 = eye(4);
        for k2 =1:4
            rho1(k2,k2) = D(k2);
        end
        rho_evol{j+1} = rho1;
        chi_2 =0;
    else
    U = U_target(:,:,j)*U;
    rho_evol{j+1} = U_target(:,:,j)*rho_evol{j}*U_target(:,:,j)';
    end
end
rho_evol{j+1}
PrettyPlotMatrix3(rho_evol{12}, 100,0)
%}
%%%%%%%%%%%%                                                                                        %%%%%%%%%%%%


% % -->> % PREPARING THE INPUT STRUCTURE PulseMod

% input preparation
    % writing variables in a structure (here this is not an issue, since the parallelization only occurs later )
    for j =1:length(QbitOrder)
        PulseMod.shape_obs{j} = make_gassian_profile(np); 
        % % -->> %  SEE that one has to define a initial pulse shape, thats because there are multiple options to 
        % % -->> % optimize below, some only find the best power, others shape, etc
        PulseMod.shape_dec{j} = make_gassian_profile(np);
        PulseMod.width_obs{j} = pulse_width;
        PulseMod.width_dec{j} = pulse_width;
        PulseMod.np_obs{j} = np; 
        PulseMod.np_dec{j} = np; 
        PulseMod.power_bound_obs{j} = power_bound_obs;
        PulseMod.power_bound_dec{j} = power_bound_dec;
        PulseMod.QbitOrder{j} = QbitOrder{j};
        PulseMod.Par = Par;
        PulseMod.U_target{j} = U_target(:,:,j);
        PulseMod.s_a_obs{j} = s_a;
        PulseMod.s_p_obs{j} = s_p;
        PulseMod.s_a_dec{j} = s_a; %s_a;
        PulseMod.s_p_dec{j} = s_p; %s_p;   
        if  QbitOrder{j}{1} == "int" % adjustments for the free evolution
            %if j <6            
            if j>3
                pw_free = 2324/2;%floor(2324*3/2);
            else
                %pw_free = 2324; % 2324*3;
                pw_free = floor(2324*(50-k)/50); % 2324*3;
            end
            np_free = pw_free;
            PulseMod.width_obs{j} = pw_free;
            PulseMod.width_dec{j} = pw_free;
            PulseMod.np_obs{j} = np_free;
            PulseMod.np_dec{j} = np_free;
            PulseMod.shape_dec{j} = make_gassian_profile(np_free);
            PulseMod.shape_obs{j} = make_gassian_profile(np_free);
        end
    end
    

% % -->> % NOW we are ready to run the optimization, there are many possible paths
% %% runnning the optimizations. 

%%%%%%%%%%%%%%%%%% 1) % % -->> % If one wants to find only the best power for a pulse in uncomment   
 % to find bestpower only:
    %[PowMatWOPhase, chi_2_wo_phase] = find_pulse_bestpowermax_PartMultiGate_Het(PulseMod, rho_evol{1}, acc, L_Spec, L_Mol);
% % -->> % Optimal power is stored in PowMat
    
 % to retrieve the generated state
    %[Ut1_opt, rhot1_opt] = get_UAndRho_Het(PulseMod, rho_evol{1},PowMatWOPhase, L_Spec,L_Mol, isfullmod);
    
%%%%%%%%%%%%%%%%%% 2) % % -->> % If one wants to find only the modulation for the phase uncomment below 
 % to find only modulated phase
    %[PowMat, PulseMod, chi_2gauss_FSPhase, rho_evol] = ...
    %    get_PowerAndFSModPhase_PartMultigate_Het(PulseMod, rhoini, verbose, L_Spec, L_Mol);
    %isfullmod = 0;
    %[Ut2_opt, rhot2_opt] = get_UAndRho_Het(PulseMod, rhoini,PowMat, L_Spec,L_Mol, isfullmod); % to retrieve the state
% % -->> % Optimal power is stored in PowMat, and PulseMod optimal pulse. 

%%%%%%%%%%%%%%%%%% 3) % % -->> % If one wants to find the full Fourier Series modulation, i.e., 
% amplitude and phase modulation
 %full modulation
    [PulseMod, chi_2, rho_evol] = get_FSMod_PartMultigate_Het(PulseMod,rho_evol{1}, L_Spec,L_Mol, verbose);
   
%% ONE CAN USE THE FOLLOWING SCRIPT TO PRINT TO .txt the pulses when parallelization is on
 % FILE PRINTING FOR PARALELIZATION

%for j =1:length(PulseMod.QbitOrder)
% printing on file the results - ready for the NMR
path = '/home/user/nmr/Matlab_codes/Modulated_Pulses/NMR_PulseShapes/NMR_PulseShapes_output/CNOT2';
which_target = 1 +2*r1; 
%Save2File([], chi_2, s_a, s_p, np, pulse_width, iter,path, shape_mod*0, shape_mod, which_target, j);
Save2File([], chi_2, s_a, s_p, np, pulse_width, iter,path, PulseMod.shape_dec{j},PulseMod.shape_obs{j}, which_target, j);
%end
%{ % Printing development
fprintf('\n Last loop, pw = %d, np = %d \n', pulse_width, np);
%end

%% Finding the new unitary
% % -->> % The following code should be used to do the "reverse engineering" i.e.,
% % -->> % have the pulse but wants to find the unitary

U_total_gen=1;
rho_evol_total = rho_evol;
free_evol = ones(pulse_width,3);
power = 1e-6;
%pw_free = 2324;
free_evol(1:pulse_width,1:2) = free_evol(1:pulse_width,1:2)*1e-4;
for j=1:length(QbitOrder)
    if j~=2
        shape_mod_obs = PulseMod.shape_obs{j};
        shape_mod_dec = PulseMod.shape_dec{j};
        [U_opt, ~] = simshapedpulse2_noPar(shape_mod_obs,shape_mod_dec,pulse_width,pulse_width, ...
                    0,0,maxpower_obs,maxpower_dec,rho_evol_total{j});% , L_Spec, L_Mol);

    else
        [U_opt, rho_opt] = simshapedpulse2_Par(free_evol, free_evol, pw_free, pw_free,0, 0,power,power,rho_evol_total{j}, L_Spec, L_Mol); 
    end
    U_target_total = U_target(:,:,j);
    chi_2 = sum(sum( real(U_opt - U_target_total).^2 + imag(U_opt - U_target_total).^2 ));
    U_total_gen = U_opt*U_total_gen;
    % % -->> %  One can calculate then the Fidelities with fidelmat
    fid = fidelmat(rho_T,rho_g);
end

isfullmod = 1;
%rhoini = kron(prj{1},prj{2});
PowMat = zeros(2,11);
[Ut_opt, rhot_opt] = get_UAndRho_Het(PulseMod, rhoini,PowMat, L_Spec,L_Mol, isfullmod);


%% Plotting pulse shapes
% % -->> %  this code snippet calls the plotPulseFFT() function with an offset
% % -->> %  This function plots the Fourier Transform of the pulse.
% % -->> %   one can use the PulseMod to get the shape or define a function

Shape = PulseMod.shape_obs{1};

t= linspace(0, length(Shape(:,3))*1e-6, length(Shape(:,3))/1.2);
Shape2 = cos((80000*t) +2000000*1e-6) +1;% + sin(50000*t +100000*1e-6) ;
%plotPulseFFT(Shape2, [], t, 1,1, mol.dq)
t= linspace(0, length(Shape(:,3))*1e-6, length(Shape(:,3)));
plotPulseFFT(Shape(:,2) ,Shape(:,1),t,1,1, mol.dq - 6000)

%% EXAMPLE 2: modulating for multiQbitgates decomposed unitaries in Homonuclear systems
% % -->> %  Same as the code above for composing Unitaries and finding gates on defined circuit
% % -->> %  But now it is for a Homonuclear system, note that all functions have an "Hom" at the 
% % -->> %  end of their names

% gate sequence and acting qbits:
trifluor1 % defining molecule
global mol spectro
[L_Spec, L_Mol] = fetch_global_val(spectro, mol);
% defining initial state
rhoini = kron(kron(prj{1},prj{1}),prj{1});


% defining Qubit Order
QbitOrder = {{"sing",0}, {"sing",1}, {"sing",2}}; %,{"int"},{"sing",2}};
%U_target(:,:,1) = [kron(expm(-1i*pi/4*sy),kron(id,id))];%[kron(id,kron(expm(-1i*pi/2*sy),expm(-1i*pi/2*sy)))];
% storting unitaries as as (U, n) 3D matrix
U_target_1 = [kron(kron(expm(-1i*pi/2*sx),id),id)];
U_target_2 = [kron(kron(id,expm(-1i*pi/4*sx)),id)];
U_target_3 = [kron(kron(id,id),expm(-1i*pi/8*sx))];
U_target(:,:,1) = U_target_1;
U_target(:,:,2) = U_target_2;
U_target(:,:,3) = U_target_3;


%defining pulse variables only one channel now = no obs/dec
pulse_width = 200; np=200;
maxpower = 20*1e3;
minpower = 1e3;
power_bound_obs = [minpower maxpower];
power_bound_dec = [minpower maxpower]; % can be declared for redundancy and avoid extra code editing

% declaring optimization variables
acc = 1e-4;
s_a = 4;
s_p = 14;
Par = [2 200]; % Defining the type of parallelization, "[0 10]" means no paralelization, "[1 10]" means 
% parallelization defined outside the function as a parfor (paralelized for loop) with 10 paralel runs, and "[2 10]" parallelization
% inside as a Multistart() MATLAB function, with 10 paralel runs.,
verbose = 1; %verbose is 0 if you don't want to print anything, and 1 to print development
power = maxpower;
pulse_width_Vec = [pulse_width pulse_width pulse_width]; 
offset = 0; % here one can prior to optimization, define an offset to frequency prior to optimization


% writing variables in a structure (here this is not an issue, since the parallelization only occurs later )
for j =1:length(QbitOrder)
    PulseMod.shape_obs{j} = make_gassian_profile(np);
    PulseMod.shape_obs{j}(:,1) = 0*ones(1,np);
    PulseMod.shape_dec{j} = make_gassian_profile(np);
    %PulseMod.shape{j}(:,1) = PulseMod.shape{j}(:,3)*90*(0.8)^(j-1);
    PulseMod.width_obs{j} = pulse_width_Vec(j);
    PulseMod.width_dec{j} = pulse_width_Vec(j);
    PulseMod.np_obs{j} = np; %; floor(np*2^(1-j));
    PulseMod.np_dec{j} = np; %; floor(np*2^(1-j));
    PulseMod.power_bound_obs{j} = power_bound_obs;
    PulseMod.power_bound_dec{j} = power_bound_dec;
    PulseMod.QbitOrder{j} = QbitOrder{j};
    PulseMod.Par = Par;
    PulseMod.U_target{j} = U_target(:,:,j);
    PulseMod.s_a_obs{j} = s_a;
    PulseMod.s_p_obs{j} = s_p;
    PulseMod.s_a_dec{j} = s_a;
    PulseMod.s_p_dec{j} = s_p;
    PulseMod.offset{j} = offset;
            if  QbitOrder{j}{1} == "int" % adjustments for the free evolution
            %if j <6            
            if j>1
                pw_free = 2324/2;%floor(2324*3/2);
            else
                pw_free = 2324; % 2324*3;
            end
            np_free = pw_free;
            PulseMod.width_obs{j} = pw_free;
            PulseMod.width_dec{j} = pw_free;
            PulseMod.np_obs{j} = np_free;
            PulseMod.np_dec{j} = np_free;
            PulseMod.shape_dec{j} = make_gassian_profile(np_free);
            PulseMod.shape_obs{j} = make_gassian_profile(np_free);
        end
end


%%%%%%%%%%%%%%%%%% % % -->> % There are 3 optimization options

%%%%%%%%%%%%%%%%%% 1) % % -->> % when we are interested only on finding the best power setting
    %[power_Vec, chi_2] = find_pulse_bestpowermax_PartMultiGate_Hom(PulseMod, rhoini, acc, L_Spec, L_Mol, verbose);

%%%%%%%%%%%%%%%%%% 2) % % -->> % Find best power for amplitude, and best Fourier Series Modulation for Phase
    % example gaussian with phase modulation
    %[power_Vec, PulseMod, chi_2_gauss_FSPhase, rho_evol] = get_PowerAndFSModPhase_PartMultigate_Hom(PulseMod, rhoini, verbose, L_Spec, L_Mol);
       
%%%%%%%%%%%%%%%%%% 3) % % -->> % Find best Fourier Series modulation for both amplitude and phase
    %full modulation
    [PulseMod , chi_2_fullFS, rho_evol] = get_FSMod_PartMultigate_Hom(PulseMod, rhoini, L_Spec, verbose);
    focusfrequency(0,0); % to reset offset to zero

% % -->> % For only one gate one can optionally use:
    %[shape_mod_obs,chi_2] = get_FSMod_1q_Hom(s_a, s_p, np, pulse_width, power, local_U_target, rho_evol, L_Spec,L_Mol, verbose, Par);



%% Plotting routines and making GIFs for dynamics vizualizations
% % -->> % save the evolved state in a rho_evol_togif cell
%PrettyPlotMatrix3(rhoini,10,0);
%PrettyPlotMatrix3(target,30,0);
    % routine for saving gifs
    clear rhotoplot1
    np=11;
    pw =201;
    dtau = pw/np;
    for k=1:4
        rhotoplot1{k} = rho_evol_togif{k};
        pwd = '/home/user/nmr/Matlab_codes/Modulated_Pulses/NMR_PulseShapes/Plots/';
        filename = append('BlochPlot_', 'processname', num2str(k), '.gif');
        PlotGIF_Rho_inTime(rhotoplot1, np, dtau, pwd, filename)
    end

%% PULSE FILE OUTPUTING AS .txt FILES

% printing on file the results - ready for the NMR
path = '/home/user/nmr/Matlab_codes/Modulated_Pulses/NMR_PulseShapes/NMR_PulseShapes_output/Clorophorm-200us_';
which_target = 0; %i1;
%shape_mod_dec = zeros(np,3);
iter = 1;
%Save2File(Coeff_opt, chi_2, s_a, s_p, np, pulse_width, iter,path, shape_mod_dec, shape_mod_obs, which_target);
Save2File([], chi_2_fullFS(1), s_a, s_p, np, pulse_width, iter,path, PulseMod.shape_dec{1},PulseMod.shape_obs{1}, which_target);
%{ % Printing development
disp('Current loop: [iter, np, U_target(i1)]')
