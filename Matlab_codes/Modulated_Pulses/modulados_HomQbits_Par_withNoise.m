%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulacao de Pulsos modulados
% Sintax Simulacao dos Pulsos e Evolucao Livre
% rgpulse(pw1,ang1,phase1,roha)
% decrgpulse(pw2,ang2,phase2,roha)
% simpulse(pw1,pw2,ang1,ang2,phase1,phase2,roha)
% pw# tempo do pulso de 90, ang# angulo de rotacao, 
% phase# fasedo pulso (0 = X, 90 = Y, 180 = -X, 270 = -Y),
% rhoa estado antes do pulso
% delay((tempo/JJ),roha)
% tempo em unidades de 1/J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SEE TUTORIAL_ModulatedPulses.m for instructions %%%%


% clearing memory
clc;
clear all;
close all;
tic
% fetching objects
spectropar;
cloroformio;

global mol spectro
% Making local variables for the pulses when parallel
[L_Spec, L_Mol] = fetch_global_val(spectro, mol);

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

% NUMERICAL CONSIDERATIONS:
% % Estimativas freq do Hamiltoniano transverso.
% pw90 = 12.7; % tempo pw90 em micro-segundos
% % H_rf = 2*pi*nu*(1/2)*sigma_y
% nu = 1/(4*pw90*1e-6); % freq. do hamiltoniano transverso em Hz para 58 dB em torno de 20 kHz
% % Fator giromagnetico do proton gp = 2.675222005(63)×10**8 rad⋅s−1⋅T−1
% giro = 2.675222*1e8;
% campo_efetivo_B1 = 2*pi*nu/giro; % da ordem de 4.6973e-04 T-1 aprox. 4.7 Gauss (1 Gauss[G] = 0.0001 Tesla [T])
% % 
% % PRIOR EXPERIMENTAL DATA:
% estimativa 27/05/2023: pw90 = 11.6375. -->  nu = 21482.0 Hz;
% pw = 11.4, p1 =10.325 --> nuw = 21930, nu1 = 24213
% estimativa 13/09/2023: 
% pw90 = 11.8875 --> nu = 21030.0 Hz , p1=10.65 -- > nu1 = 23474.0




% ------------------------------------------------------------------------ %
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< %
% ------------------------------------------------------------------------ %

%%%%%%%%%%%%%%% starting functions %%%%%%%%%%%

function [L_spectro, L_mol] = fetch_global_val(spectro, mol) 
    % saving locally global variables to get around paralelization restrictions
    % Please notice how indexes are related to structure methods
    %nspin = mol.nspin();
    L_spectro(1) = spectro.rof1;
    L_spectro(2) = spectro.rof2;
    L_mol(:,:,1) = mol.Hint;
    L_mol(:,:,2) = mol.Hzee;
    for j =1:mol.nspin()
        L_mol(:,:,2*j+1) = mol.Ix{j};
        L_mol(:,:,2*j+2) = mol.Iy{j};
    end

end


%%%%%%%%%%%%%%% plotting / file handling functions %%%%%%%%%%%


function [phase, amplitude, timesteps] = fetch_shape_from_file(path)
    fileID = fopen(path, 'r');
    % Skip commented lines
    line = fgets(fileID);
    while line(1) == '#'
        line = fgets(fileID);
    end
    % Read the data
    data = textscan(fileID, '%f %f %f');
    % Close the file
    fclose(fileID);
    % Extract phase, amplitude, and timesteps
    phase = data{1};
    amplitude = data{2};
    timesteps = data{3}; % You may or may not use this, depending on your needs
end

function PlotGIF_Rho_inTime(rho_evol, Ntimesteps, timestep_size, path, filename)
    % Makes a gif from the plots of PrettyPlotMatrix3
    % rho_evol = should be a Cell array containing the matrices rho_evol(timestep), 
    % Ntimesteps = is number of timesteps; timestep_size = number of points could be smaller than timesteps 
            % (ex: one rho_evol for each 10 seconds, so that time flows t --> t + timestepsize, for each 0 < j < Ntimesteps
    % path should be = '/path_to_saving_directory'; filename = string name of the file
    filename = append(path, '/' , filename);
    np = Ntimesteps;
    dtau = timestep_size; %this just prints the timestep, it doesn't simulate it 

    for j = 1:np %np % loop for saving plots
        %PrettyPlotMatrix3(rho_evol{j},20,0);
        for jj =1:j
            rho_evol_subset{jj} = rho_evol{jj};
        end
        PlotBloch(rho_evol_subset, 20,30,[0.25,0,0.25],false);

        annotation("textbox",[.07 .6 .3 .3], "String", append("time (us) = ", int2str(j*dtau)), "FitBoxToText", "on" );
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
            if j == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',0.8);
            elseif j>1
                imwrite(imind,cm,filename,"gif","WriteMode","append","DelayTime",0.8);
                if j==np
                    imwrite(imind,cm,filename,"gif","WriteMode","append","DelayTime",2.2);
                    imwrite(imind,cm,filename,"gif","WriteMode","append","DelayTime",2.2);
                end
            end
        delete(findall(gcf,'type','annotation'))
    end
end

function res = plot_AmpAndPhase(nfig,shape_mod, np, pulse_width)
    %plots routines for visualization
    % routine for x axis: remember shape_mod(:,3) = timesteps
    shape_mod(:,3) = shape_mod(:,3)*(pulse_width/np); % set 3rd column = timesteps
    for j =2:length(shape_mod(:,3)) %iterate on timesteps cummulatively
        shape_mod(j,3) = shape_mod(j,3) + shape_mod(j-1,3);
    end

    figure(nfig);
    subplot(1,2,1)
    %plotting amp
    plot(shape_mod(:,3),shape_mod(:,2))
    ylim([0 1023]);
    title('Fourier series Pulse Amplitude')
    subplot(1,2,2) 
    %OBS: this plot removes "discontinuity" points, because
    %matlab doesn't understands mod() plots
    Phase2Plot = shape_mod(:, 1);
    Disc2Nan = abs(diff(shape_mod(:, 1)))>300;
    for j = 1:size(Disc2Nan)
        if Disc2Nan(j) == 1
            Phase2Plot(j+1) = NaN;
        end
    end
    plot(shape_mod(:,3), Phase2Plot,  "LineWidth", 2)
    %xticklabels({'0','1.5', '3.0'})
    ylim([0,360])
    title('Fourier series Pulse Phase')
end

function SavingShape = Save2File(Coeff,chi_2, s_a, s_p, np, pulse_width, iter, path, shape_mod_dec, shape_mod_obs, target, step)
    %function for saving the pulse shapes in matlab format. 
    % The saved archive should look like:
    %   %%%% commments defining constants
    %   phase_column     amplitude_column     time_column
    % OBS: Check machine precision before writing to file, otherwise final resolution won't correspond to file print
    
    if target == 1
        str_target = 'Its a Rotx(90) Id gate';
        gate = "(Rotx(90)Xid)_";
    elseif target==2
        str_target = 'Its a 1 tensorprod sigma_y gate';
        gate = "(1Xsy)_";
    elseif target == 3
        str_target = 'Its a Roty(90) Id gate';
        gate = "(Roty(90)Xid)_";
    elseif target == 0
        str_target = 'Its a sigma_x gate tensorprod 1 tensorprod 1';
        gate = "(sxX1X1)_";
    
    end
    
    % for machine precision:

    % Turning numbers to strings
    iter = int2str(iter); np = "NP-" + int2str(np) + "_"; 
    S = " s_a,s_p = " + int2str(s_a) + "," + int2str(s_p);
    Total_time = append('Pulse Width = ',num2str(pulse_width));
    chi_2 = "Chi_2 = " + num2str(chi_2);
    numCoef = "SA" + int2str(s_a) + "SP" + num2str(s_p) + "_";
    step = append('step', num2str(step), "_");
    
    target_header = append('# Gate tipe \n# ' , str_target, '\n# \n');    
    header = append('# Some important values: \r\n# ', chi_2,' , \n#', S , ' , ', Total_time, '\r\n' ,'# Coefficients : \n');
    header2 = append('# \n# Phase Amplitude Timesteps  \n# \n');    
    
    %commenting Dec out
    
    name= append(path,'Shape_', gate, step,'PW-',num2str(pulse_width), np, numCoef,'_Opt', iter,'_Dec', '.RF');
    fid = fopen(name, 'w');
    fprintf(fid, target_header);
    fprintf(fid, header);
    fprintf(fid, '# %f %f %f \r\n', [Coeff]');
    fprintf(fid, header2);
    fprintf(fid, '%f     %f     %f \r\n', [shape_mod_dec]');
    fid = fclose(fid);
    
    name= append(path,'Shape_', gate,step,'PW-',num2str(pulse_width),np, numCoef,'_Opt', iter, '_Obs', '.RF');
    fid = fopen(name, 'w');
    fprintf(fid, target_header);
    fprintf(fid, header);
    fprintf(fid, '# %f  %f  %f \r\n', [Coeff]');
    fprintf(fid, header2);
    fprintf(fid, '%f     %f     %f \r\n', [shape_mod_obs]');
    fid = fclose(fid);

    
end

%%%%%%%%%%%%% pulse functions %%%%%%%%%%%

function Coeff_ini = Coeff_Guess(s_a,s_p,flag)
    % generates a matrix Coeff_ini such that it is a concatenation of the coefficients that define an initial
    % Fourier series for amplitude and phase of a given pulse
    % Nq int = 1 % number of qubits in gate
    % also: size(flag) = Nq
    arguments
        s_a double = 1 %amplitude series rank
        s_p double = 1 %phase_series rank
        flag (1,:) = 'Rweighted' % options are 'weighted','ones', 'zeros', 'singlemode', 'random', 'Rweighted'
    end
    
    Nq = size(flag,2);
    Coeff_ini=[];
    for q = 1:Nq
        A_Coef = []; P_Coef = [];
        if flag(1,q) == "weighted"
            for j = 1:s_a
                A_Coef(j,1) = 0.5/j;
                A_Coef(j,2) = 0.07/j;
                A_Coef(j,3) = 0.0;
            end
            for j = 1:s_p
                P_Coef(j,1) = 90/j;
                P_Coef(j,2) = 0.20/j;
                P_Coef(j,3) = 0.0;
            end
        elseif flag(1,q) == "Rweighted"
            for j = 1:s_a
                A_Coef(j,1) = (1/8 + (1-1/8)*rand(1))*120/sqrt(j);
                A_Coef(j,2) = (1/8 + (1-1/8)*rand(1))*0.12;%*0.12;
                A_Coef(j,3) = (1/8 + (1-1/8)*rand(1))*36;
            end
            for j = 1:s_p
                P_Coef(j,1) = (1/8 + (1-1/8)*rand(1))*120/sqrt(j);
                P_Coef(j,2) = (1/8 + (1-1/8)*rand(1))*0.14;%0.18
                P_Coef(j,3) = (1/8 + (1-1/8)*rand(1))*36;
            end
        elseif flag(1,q) == "Rweighted2"
            for j = 1:s_a
                A_Coef(j,1) = (1/8 + (1-1/8)*rand(1))*120;%/sqrt(j);
                A_Coef(j,2) = (1/8 + (1-1/8)*rand(1))*0.04;%*0.04;
                A_Coef(j,3) = (1/8 + (1-1/8)*rand(1))*36;
            end
            for j = 1:s_p
                P_Coef(j,1) = (1/8 + (1-1/8)*rand(1))*120;%/sqrt(j);
                P_Coef(j,2) = (1/8 + (1-1/8)*rand(1))*0.08;%0.08
                P_Coef(j,3) = (1/8 + (1-1/8)*rand(1))*36;
            end
        elseif flag(1,q) == "ones"
             A_Coef = 0.5*ones(s_a,3);
             P_Coef = 0.5*ones(s_p,3);
        elseif flag(1,q) =="zeros"
             A_Coef = zeros(s_a,3);
             P_Coef = zeros(s_p,3);
        elseif flag(1,q) == "singlemode"
             A_Coef(1,:) = [0.5 0.05 pi/2];
             P_Coef(1,:) = [90/2 0.1 pi/2];
        elseif flag(1,q) == "random"
             A_Coef = rand(s_a,3);
             P_Coef = rand(s_p,3);
        end
        C_temp = cat(1,A_Coef,P_Coef);
        Coeff_ini = cat(1, Coeff_ini, C_temp); % Initial guess for Amp and phase coefs. concatenation
    end
end

function [power_Vec, PulseMod, chi_2, rho_evol_opt] = ...
             get_PowerAndFSModPhase_PartMultigate_Hom(PulseMod, rhoini, verbose, L_Spec, L_Mol);
%finds the best power, and phase modulation on each unitary applied on a circuit, assumes that the Pulse
% shapes and properties are all in PulseMod structure and that unitaries are partitioned 
% and ordered in PulseMod.U_target. 
arguments
    PulseMod struct %structure with pulse properties
    rhoini (:,:) %inicial state
    verbose double % 0 for False, 1 for true
    L_Spec %spectropar as local variable
    L_Mol (:,:,:) % Molecules Ix's and Iy's stored as a local variable
end
    rho_evol{1} = rhoini;
    rho_evol_opt{1} = rhoini;
    acc = 1e-4;

    for j =1:length(PulseMod.QbitOrder)
        if verbose ==1
            fprintf('\noptimizing gate %d of %d', j , length(PulseMod.QbitOrder))
        end
        % unpacking PulseMod
        currentQbit = PulseMod.QbitOrder{j};
        U_target = PulseMod.U_target{j};
        pw_obs = PulseMod.width_obs{j};
        power_bound_obs = PulseMod.power_bound_obs{j};
        AmpShape = PulseMod.shape_obs{j};
        s_p = PulseMod.s_p_obs{j};
        np = PulseMod.np_obs{j};
        Par = PulseMod.Par;

        switch currentQbit{1}
            case {"int", 'grad'} % free evolution or gradient
                pw = pw_obs;
                free_evol = ones(pw,3);
                free_evol(1:pw,1:2) = free_evol(1:pw,1:2)*1e-4;
                power = 1e-6;
                if currentQbit{1} == "int" %free evolution
                    %shapedpulse_Hom_Par(shape_mod,pulse_width,shape_mod(:,1),0,power,rho_ini, L_Spec, L_Mol); 
                        [U_gen, rho_opt] = shapedpulse_Hom_Par(free_evol, pw, [],0,power,rho_evol{j}, L_Spec, L_Mol); 
                        chi_2 = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 )) + abs(angle(U_gen(1,1)) - angle(U_target(1,1)));
                        if verbose ==1
                            fprintf("\nModulaton metrics (free evolution): \n chi_2 = %d , for pw = %d \n", chi_2, pw)
                        end
                        rho_evol{j+1} = U_target*rho_evol{j}*U_target';
                else % gradient                    
                    if verbose == 1
                        fprintf('running gradient on gate number %d', j)
                    end
                    D = diag(rho_evol{j});
                    rho1 = eye(4);
                    for k2 =1:4
                        rho1(k2,k2) = D(k2);
                    end
                    rho_evol{j+1} = rho1;
                    chi_2 =0;
                end
                shape_mod_obs = free_evol;
            case{'sing'}         
                % setting resonance frequency
                target_qbit = currentQbit{2};
                offset = PulseMod.offset{j};
                L_Mol = focusfrequency(target_qbit, offset);      
            
                % gaussian modulation finds phase and power
                [power_Vec(j), shape_mod_obs, chi_2(j)] =  ...
                    get_PowerAndFSModPhase_1q_Hom(AmpShape, s_p, np, pw_obs, power_bound_obs,...
                                    rho_evol{j}, U_target, acc, verbose, Par, L_Spec, L_Mol);
                  
                % getting the final [U, rhofin]
                %[U_gen, rhofin] = shapedpulse_Hom_Par(shape,pw,[],0,(maxpower_loop -r*acc*1e4),rho_ini, L_Spec, L_Mol); 
                [U_opt, ~] = shapedpulse_Hom_Par(shape_mod_obs,pw_obs,[],0,power_Vec(j),rho_evol{j}, L_Spec, L_Mol);
                %chi_22(j) = sum(sum( real(U_opt - U_target).^2 + imag(U_opt - U_target).^2 ));
                rho_evol{j+1} = U_target*rho_evol{j}*U_target';
                rho_evol_opt{j+1} = U_opt*rho_evol_opt{j}*U_opt';
                PulseMod.shape_obs{j} = shape_mod_obs;
        end
    end
end

function [power, shape_mod_obs, chi_2] = get_PowerAndFSModPhase_1q_Hom(AmpShape, s_p, np, pulse_width, power_bound, ...
                                        rhoini, U_target, acc, verbose, Par, L_Spec, L_Mol)
    
    %>>>>>>>>>>>> OBS <<<<<<<<<<<<<
    % see later if adding some power after f
    % inding the best power doesn't give better chi2's ...
    % spoiler: it does! 
    %>>>>>>>>>>>> <<<<<<<<<<<<<<<<<

    %finding best amplitude/power for given pulse width and maxpower and modulates phase to better approximate a give unitary
    %arguments below are given in a specific case of 2 qbits, but it runs for more than 2
    %paralellization runs as follows: isPar = 0 for non-parallel, 1 for parallelization defined outside the function as a parfor
    % 2 for parallelization inside as a multistart, it then has a second argument as Parruns = number of starting points
    arguments
        AmpShape (:,:) = [] % if is an .RF file, load it outside of the function
        s_p double = 8
        np double = 200
        pulse_width double = 200
        power_bound = [0, 20*1e3]
        rhoini (:,:) = [1 0; 0 0]
        U_target (:,:)  =  [kron(expm(-1i*pi/2*sx),id)];
        acc double = 1*e-4
        verbose = 1 %set to 1 if verbose, set to 0 if not verbose
        Par (1,:) = [0 44]
        L_Spec (:,:,:) = [] % spetropar local variable Only used when parallelizing
        L_Mol (:,:,:) = [] % mol local variable
    end
    isPar = Par(1);
    Par_runs = Par(2);

    %function [power, chi_2_opt] = set_180pulse_powermax_Hom(shape, pw, maxpower, rho_ini, U_target, targetQbit,acc) 
    [power, chi_2_wo_phase] = find_pulse_bestpowermax_Hom(AmpShape, pulse_width, power_bound, rhoini, U_target, acc, L_Spec,L_Mol);
    %power = power +2000;
    power_bound = [power_bound(1), power];
    
    Maxeval = (7+s_p)*floor(1/acc);
    if verbose == 1
        fprintf('\n Found best power in Hz: power = %d \n chi_2 prior to optimization: chi_2 = %d, \n searching phase ...', power, chi_2_wo_phase)
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval, ...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-5);
    else
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval,'TolX', 1e-5);
    end

    flag = "Rweighted"; % select a inicial condition 'weighted','random','zeros','ones'
    Coeff_ini = Coeff_Guess(0,s_p, flag); %number of qbits Nq = size(flag,1)
    %Coeff_ini(:,2)/50;
    
    if isPar== 2
        % some parallelization configurations
        ms = MultiStart('FunctionTolerance',2e-15,'MaxTime', 3600,'UseParallel',true);
        problem = createOptimProblem('fmincon','x0',Coeff_ini,... 
                'objective', @(Coeff) find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power,pulse_width, rhoini, L_Spec, L_Mol,isPar-1), ...
               'lb',-12/5*Coeff_ini,'ub',12/5*Coeff_ini);  %8/5Coeff_ini
        %runs the ms problem with Par_runs cpus
        Coeff_opt = run(ms,problem,Par_runs);
    elseif isPar == 1
        %function find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rho_ini, L_Spec, L_Mol, isPar)
        %[Coeff_opt,chi_2,exitflag,output] = fminsearch(@(Coeff) ...    
        %    find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar), ...
        %    Coeff_ini,options) ;% optimization with fminsearch
        [Coeff_opt,chi_2,exitflag,output] = fminsearchbnd(@(Coeff) ...
            find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar), ...
            Coeff_ini,-12/5*Coeff_ini,12/5*Coeff_ini,options) ;% optimization with fminsearch with a taylored function, doesn`t work very well with upper bounds...
    elseif isPar ==0
        [Coeff_opt,chi_2,exitflag,output] = fminsearchbnd(@(Coeff) ...
            find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar), ...
            Coeff_ini,-12/5*Coeff_ini,12/5*Coeff_ini,options) ;% optimization with fminsearch with a taylored function, doesn`t work very well with upper bounds...
    end
    shape_mod_obs = make_FourierProfile([], Coeff_opt, np);
    shape_mod_obs(:,2) = AmpShape(:,2);
    if isPar ==2
        chi_2 = find_phase_FSModulation_loss1q_Hom(U_target,Coeff_opt,AmpShape,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar-1);
    else
        chi_2 = find_phase_FSModulation_loss1q_Hom(U_target,Coeff_opt,AmpShape,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar);
    end
    if verbose == 1
        fprintf('\n found best modulation with chi^2 = %d', chi_2)        
        plot_AmpAndPhase(30, shape_mod_obs,np, pulse_width)
    end
end

function [power_Vec, chi_2] = find_pulse_bestpowermax_PartMultiGate_Hom(PulseMod, rhoini, acc, L_Spec,L_Mol,verbose)
% finds the best power on each unitary applied on a circuit, assumes that the Pulse shapes and properties are all in PulseMod
% structure and that unitaries are partitioned and ordered in PulseMod.U_target. 
    arguments
        PulseMod struct %Pulse properties
        rhoini (:,:) %inicial state
        acc double % precision on power, usual 1e-4
        L_Spec (:,:,:) % spetropar as local variable
        L_Mol (:,:,:) % molecules Ix's, Iy's, Hzee and Hint stored locally
        verbose = 0
    end
    rho_evol{1} = rhoini;
    rho_evol_opt{1} = rhoini;
    acc = 1e-4;

    for j =1:length(PulseMod.QbitOrder)
        if verbose ==1
            fprintf('\noptimizing gate %d of %d', j , length(PulseMod.QbitOrder))
        end
        % unpacking PulseMod
        currentQbit = PulseMod.QbitOrder{j};
        U_target = PulseMod.U_target{j};
        pw_obs = PulseMod.width_obs{j};
        power_bound_obs = PulseMod.power_bound_obs{j};
        AmpShape = PulseMod.shape_obs{j};
        np = PulseMod.np_obs{j};
        Par = PulseMod.Par;

        switch currentQbit{1}
            case {"int", 'grad'} % free evolution or gradient
                pw = pw_obs;
                free_evol = ones(pw,3);
                free_evol(1:pw,1:2) = free_evol(1:pw,1:2)*1e-4;
                power = 1e-6;
                if currentQbit{1} == "int" %free evolution
                    %shapedpulse_Hom_Par(shape_mod,pulse_width,shape_mod(:,1),0,power,rho_ini, L_Spec, L_Mol); 
                        [U_gen, rho_opt] = shapedpulse_Hom_Par(free_evol, pw, [],0,power,rho_evol{j}, L_Spec, L_Mol); 
                        chi_2 = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 )) + abs(angle(U_gen(1,1)) - angle(U_target(1,1)));
                        if verbose ==1
                            fprintf("\nModulaton metrics (free evolution): \n chi_2 = %d , for pw = %d \n", chi_2, pw)
                        end
                        rho_evol{j+1} = U_target*rho_evol{j}*U_target';

                else % gradient                    
                    if verbose == 1
                        fprintf('running gradient on gate number %d', j)
                    end
                    D = diag(rho_evol{j});
                    rho1 = eye(4);
                    for k2 =1:4
                        rho1(k2,k2) = D(k2);
                    end
                    rho_evol{j+1} = rho1;
                    chi_2 =0;
                end
                power_Vec(j) = 0;
                shape_mod_obs = free_evol;
            case{'sing'}             

                % setting resonance frequency
                target_qbit = currentQbit{2};
                offset = PulseMod.offset{j};
                L_Mol = focusfrequency(target_qbit, offset);    

                %only pulse no modulation
                
                [power, chi_2_wo_phase] = find_pulse_bestpowermax_Hom(AmpShape, pw_obs, power_bound_obs, rho_evol{j}, U_target, acc, L_Spec,L_Mol);
                power_Vec(j) = power;
            
                % getting the final [U, rhofin]
                %[U_gen, rhofin] = shapedpulse_Hom_Par(shape,pw,[],0,(maxpower_loop -r*acc*1e4),rho_ini, L_Spec, L_Mol); 
                [U_opt, ~] = shapedpulse_Hom_Par(AmpShape,pw_obs,[],0,power,rho_evol{j}, L_Spec, L_Mol);
                chi_2(j) = sum(sum( real(U_opt - U_target).^2 + imag(U_opt - U_target).^2 ));
                rho_evol{j+1} = U_target*rho_evol{j}*U_target';
        end    

    end
    
end

function [power, chi_2_opt] = find_pulse_bestpowermax_Hom(shape, pw, power_bound, rho_ini, U_target, acc,L_Spec, L_Mol) 
    % from a pulse_width especification and a fixed regular shape modulation (e.g. hard.RF, gauss.RF, isech.RF)
    % finds the best amplitude to do a pi pulse.
    % this code doesn't run the actual maxpower, since it is assumed that it is experimentally forbiden
    % default should be for 180 degrees on first qbit
    arguments 
        shape (:,3) % loaded from an '.RF' file
        pw double
        power_bound (1,2) double = [0 2*1e4] %works for 2*1eN
        rho_ini (:,:) double = eye(8)/8 % should be density matrix, usually in state |0^N>
        U_target (:,:) double = kron([0 1; 1 0], kron(eye(2), eye(2)))
        acc double = 1e-4 %must be powers of 10
        L_Spec (:,:,:) = [] % spetropar local variable Only used when parallelizing
        L_Mol (:,:,:) = [] % mol local variable
    end
    minpower = power_bound(1,1);
    maxpower = power_bound(1,2);
    % routine for non rounded maxpower*1e-N
    lexp = floor(log10(maxpower));
    hexp = floor(maxpower*10^(-lexp+1));
    hexp = hexp/10;
      
    %setting initial values
    Npulse = floor(1/acc);
    power = maxpower;
    chi_2_opt = 100;


    if Npulse>10
        [power, chi_2_opt] = find_pulse_bestpowermax_Hom(shape, pw, power_bound,rho_ini, U_target,acc*10,L_Spec,L_Mol); 
        Nlim = 2*10 -1;
        maxpower_loop = power;
        if (power)<=maxpower*(1-acc/2) %some conditions to not exceed maximum power bound.
            maxpower_loop = maxpower_loop +maxpower*acc*10/2;
        end
        for r=1:Nlim %garantee to run +-10% maxpower over last found power
            %[U rohfin] = shapedpulse2(shape,pw,phase,ramp,power,roha)
            [U_gen, rhofin] = shapedpulse_Hom_Par(shape,pw,[],0,(maxpower_loop -r*acc*1e4),rho_ini, L_Spec, L_Mol); 
            chi_2_r = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 ));
            %chi_2_r = sum(sum(abs(abs(U_gen) - abs(U_target))));
            if chi_2_r<chi_2_opt
                p = maxpower_loop -r*acc*1e4;
                if p >= minpower
                    power = p;
                    chi_2_opt = chi_2_r;
                    if chi_2_opt < acc
                        return
                    end
                end
            end
        end
    else
        for r=1:hexp*Npulse 
            r=r/2;
            %[U rohfin] = shapedpulse2(shape,pw,phase,ramp,power,roha)
            [U_gen, rhofin] = shapedpulse_Hom_Par(shape,pw,[],0,(1-r*acc)*maxpower,rho_ini, L_Spec, L_Mol); 
           %chi_2_r = sum(sum(abs(abs(U_gen) - abs(U_target))));
            chi_2_r = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 ));
            if chi_2_r<chi_2_opt
                p = (1-r*acc)*maxpower;
                if p >= minpower
                    power = p;
                    chi_2_opt = chi_2_r;
                    if chi_2_opt < acc
                        return
                    end
                end
            end
        end
    end
end

function [PowMat, PulseMod, chi_2gauss_FSPhase, rho_evol] = ...
    get_PowerAndFSModPhase_PartMultigate_Het(PulseMod, rhoini, verbose, L_Spec, L_Mol)
    
        % ---------------------% DOESNT SIMULATE "grad" for now...
    rho_evol{1} = rhoini;
    rho_evol_opt{1} = rhoini;
    acc = 1e-4;

    for j =1:length(PulseMod.QbitOrder)
        % unpacking PulseMod
        currentQbit = PulseMod.QbitOrder{j}{1};
        U_target = PulseMod.U_target{j};
        AmpShape_obs = PulseMod.shape_obs{j};
        AmpShape_dec = PulseMod.shape_dec{j};
        pw_obs = PulseMod.width_obs{j};
        pw_dec = PulseMod.width_dec{j};
        power_bound_obs = PulseMod.power_bound_obs{j};
        power_bound_dec = PulseMod.power_bound_dec{j};
        s_pVec(1) = PulseMod.s_p_obs{j};
        s_pVec(2) = PulseMod.s_p_dec{j};
        np = PulseMod.np_obs{j};
        Par = PulseMod.Par;
    

        if currentQbit =="int" % free evolution
            pw = pw_obs;
            free_evol = ones(pw,3);
            free_evol(1:pw,1:2) = free_evol(1:pw,1:2)*1e-4;
            power = 1e-6;
            [U_gen, rho_opt] = simshapedpulse2_Par(free_evol, free_evol, pw, pw,0, 0,power,power,rho_evol{j}, L_Spec, L_Mol); 
            chi_2 = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 ));
            shape_mod_obs = free_evol;
            shape_mod_dec = free_evol;
            power_obs = power;
            power_dec = power;
            if verbose ==1
                fprintf("\nModulaton metrics (free evolution): \n chi_2 = %d , for pw = %d \n", chi_2, pw)
            end
        elseif currentQbit == "sim"
        [power_obs, power_dec, shape_mod_obs, shape_mod_dec, chi_2] = ...
            get_PowerAndFSModPhase_2q_Het(AmpShape_obs, AmpShape_dec,s_pVec, np, [pw_obs pw_dec], ...
            power_bound_obs,power_bound_dec, rho_evol{j}, U_target, acc,verbose,Par, L_Spec, L_Mol);
        elseif currentQbit =="sing"
            chn = PulseMod.QbitOrder{j}{2};
            [power, shape_mod, chi_2] = ...
            get_PowerAndFSModPhase_1q_Het(AmpShape_obs, AmpShape_dec,s_pVec(chn+1), np, [pw_obs pw_dec], ...
            power_bound_obs,power_bound_dec, rho_evol{j}, U_target, acc,verbose,Par, L_Spec, L_Mol,chn);
            if chn ==0
                power_obs = power;
                power_dec = 1e-6;
                shape_mod_obs = shape_mod;
                shape_mod_dec = shape_mod*1e-6;
            else
                power_dec =power;
                power_obs = 1e-6;
                shape_mod_dec = shape_mod;
                shape_mod_obs = shape_mod*1e-6;
            end
        end
        rho_evol{j+1} = U_target*rho_evol{j}*U_target';
        
        chi_2gauss_FSPhase(j) = chi_2;
        PowMat(1,j) = power_obs;
        PowMat(2,j) = power_dec;
        PulseMod.shape_obs{j} = shape_mod_obs;
        PulseMod.shape_dec{j} = shape_mod_dec;

    end
end

            %[power, shape_mod, chi_2] = ...
            %get_PowerAndFSModPhase_1q_Het(AmpShape_obs, AmpShape_dec,s_pVec(chn+1), np, [pw_obs pw_dec], ...
            %power_bound_obs,power_bound_dec, rho_evol{j}, U_target, acc,verbose,Par, L_Spec, L_Mol,chn);

function [power, shape_mod, chi_2] = get_PowerAndFSModPhase_1q_Het( ...
    AmpShape_obs, AmpShape_dec, s_p, np, ...
    pwVec, power_bound_obs, power_bound_dec, rhoini, U_target, acc, verbose, Par, L_Spec, L_Mol, chn)
    %finding best amplitude/power for given pulse width and maxpower and modulates phase to better approximate a give unitary
    %arguments are given in a specific case of 2 qbits
    arguments
        AmpShape_obs (:,:) = [] % if is an .RF file, load it outside of the function
        AmpShape_dec (:,:) = [] % if is an .RF file, load it outside of the function
        s_p double = 1 % vector form for [s_p_obs s_p_dec]
        np double = 200
        pwVec (1,:) = [200 200] % vector form for [pulse_width_obs pulse_width_dec] 
        power_bound_obs = [0, 20*1e3]
        power_bound_dec = [0, 20*1e3]
        rhoini (:,:) = [1 0; 0 0]
        U_target (:,:)  =  [kron(expm(-1i*pi/2*sx),id)];
        acc double = 1*e-4
        verbose = 1 %set to 1 if verbose, set to 0 if not verbose
        Par (1,:) = [0 44]
        L_Spec (:,:,:) = [] % spetropar local variable Only used when parallelizing
        L_Mol (:,:,:) = [] % mol local variable
        chn double = 0 %0 or 1
    end
    % getting variables
    isPar = Par(1);
    Par_runs = Par(2);
    pw_obs = pwVec(1);
    pw_dec = pwVec(2);
    power_dec = power_bound_dec(1,2);
    power_obs = power_bound_obs(1,2);

    %finding powers with constant phase
    if chn == 0;
        [power, chi_2] = find_pulse_bestpowermax_Het(chn, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                            power_bound_dec, rhoini, U_target, acc);
        power_obs = power;
    elseif chn == 1;
    [power, chi_2] = find_pulse_bestpowermax_Het(chn, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                            power_bound_dec, rhoini, U_target, acc);
       power_dec = power;
    end
    
    % hardcoding extra bulge in power for finding a better modulation
    % if power is already too high (above 19K), consider doing a longer pulse_width (pw_obs/dec)
    %{
    if power_obs< 19*1e3
        power_obs = power_obs+500;% +500;
    end
    if power_dec < 19*1e3
        power_dec = power_dec+500;% +500;
    end
    %}
    Maxeval = (7+(s_p))*floor(1/acc);
    if verbose == 1
        fprintf(['\n Found best power in Hz: power = %d ...' ...
            ' \n chi_2 prior to optimization: chi_2 = %d, \n searching phase ...'], power, chi_2)
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval, ...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-5);
    else
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval,'TolX', 1e-5);
    end

    flag = "Rweighted2"; % select a inicial condition 'weighted','random','zeros','ones'
    Coeff_ini = Coeff_Guess(0,s_p, flag); %number of qbits Nq = size(flag,1) %this step may be better with separated s_p's
    
    if isPar== 2
     %result = find_phase_FSModulation_loss1q_Het(chn, U_target,Coeff,AmpShape_obs, AmpShape_dec,s_p,np, ...
     %                           power_obs, power_dec, pw_obs, pw_dec, rho_ini, L_Spec, L_Mol, isPar)  
        % some parallelization configurations
        ms = MultiStart('FunctionTolerance',2e-15,'MaxTime', 3600,'UseParallel',true);
        problem = createOptimProblem('fmincon','x0',Coeff_ini,... 
                'objective', @(Coeff) find_phase_FSModulation_loss1q_Het(chn ,U_target,Coeff,AmpShape_obs, AmpShape_dec, ...
                                    s_p, np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar-1), ...
               'lb',-12/5*Coeff_ini,'ub',12/5*Coeff_ini);  %8/5Coeff_ini
        %runs the ms problem with Par_runs cpus
        Coeff_opt = run(ms,problem,Par_runs);
    elseif isPar ==1
        %function find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rho_ini, L_Spec, L_Mol)  
        %[Coeff_opt,chi_2,exitflag,output] = fminsearch(@(Coeff) ...    
        %    find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar), ...
        %    Coeff_ini,options) ;% optimization with fminsearch
        [Coeff_opt,chi_2,exitflag,output] = fminsearchbnd(@(Coeff) ...
            find_phase_FSModulation_loss1q_Het(chn ,U_target,Coeff,AmpShape_obs, AmpShape_dec, ...
                                    s_p, np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar), ...
            Coeff_ini,-12/5*Coeff_ini,12/5*Coeff_ini,options) ;% optimization with fminsearch with a taylored function, doesn`t work very well with upper bounds...
    else %no parallelization
        [Coeff_opt,chi_2,exitflag,output] = fminsearchbnd(@(Coeff) ...
            find_phase_FSModulation_loss1q_Het(chn ,U_target,Coeff,AmpShape_obs, AmpShape_dec, ...
                                    s_p, np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar), ...
            Coeff_ini,-12/5*Coeff_ini,12/5*Coeff_ini,options) ;% optimization with fminsearch with a taylored function, doesn`t work very well with upper bounds...
    end
    shape_mod = make_FourierProfile([], Coeff_opt(1:s_p, :), np);
    if chn==0
    shape_mod(:,2) = AmpShape_obs(:,2);
    chi_2 = find_phase_FSModulation_loss1q_Het(chn ,U_target,Coeff_opt,shape_mod, AmpShape_dec, ...
                                    s_p, np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar);
    elseif chn==1
    shape_mod(:,2) = AmpShape_dec(:,2);
    chi_2 = find_phase_FSModulation_loss1q_Het(chn ,U_target,Coeff_opt,AmpShape_obs, shape_mod, ...
                                    s_p, np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar);
    end
    
    if verbose == 1
        fprintf('\n found best modulation with chi^2 = %d', chi_2)        
        plot_AmpAndPhase(3, shape_mod,np, pw_obs)
    end
end

function [power_obs, power_dec, shape_mod_obs, shape_mod_dec, chi_2] = ...
        get_PowerAndFSModPhase_2q_Het(AmpShape_obs, AmpShape_dec, s_pVec, np, ...
                                    pwVec, power_bound_obs, power_bound_dec, rhoini, U_target, acc, verbose, Par, L_Spec, L_Mol)
    %finding best amplitude/power for given pulse width and maxpower and modulates phase to better approximate a give unitary
    %arguments are given in a specific case of 2 qbits
    arguments
        AmpShape_obs (:,:) = [] % if is an .RF file, load it outside of the function
        AmpShape_dec (:,:) = [] % if is an .RF file, load it outside of the function
        s_pVec (1,:) = [1 1] % vector form for [s_p_obs s_p_dec]
        np double = 200
        pwVec (1,:) = [200 200] % vector form for [pulse_width_obs pulse_width_dec] 
        power_bound_obs = [0, 20*1e3]
        power_bound_dec = [0, 20*1e3]
        rhoini (:,:) = [1 0; 0 0]
        U_target (:,:)  =  [kron(expm(-1i*pi/2*sx),id)];
        acc double = 1*e-4
        verbose = 1 %set to 1 if verbose, set to 0 if not verbose
        Par (1,:) = [0 44]
        L_Spec (:,:,:) = [] % spetropar local variable Only used when parallelizing
        L_Mol (:,:,:) = [] % mol local variable
    end
    % getting variables
    isPar = Par(1);
    Par_runs = Par(2);
    s_p_obs = s_pVec(1);
    s_p_dec = s_pVec(2);
    pw_obs = pwVec(1);
    pw_dec = pwVec(2);

    %finding powers with constant phase
    channel = 0;
    [power_obs, chi_2_wo_phase_obs] = find_pulse_bestpowermax_Het(channel, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                            power_bound_dec, rhoini, U_target, acc);
    power_bound_obs(1,2) = power_obs;
    channel = 1;
    [power_dec, chi_2_wo_phase_dec] = find_pulse_bestpowermax_Het(channel, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                            power_bound_dec, rhoini, U_target, acc);
    
    % hardcoding extra bulge in power for finding a better modulation
    % if power is already too high (above 19K), consider doing a longer pulse_width (pw_obs/dec)
    %{
    if power_obs< 19*1e3
        power_obs = power_obs+500;% +500;
    end
    if power_dec < 19*1e3
        power_dec = power_dec+500;% +500;
    end
    %}
    Maxeval = (7+(s_p_obs+s_p_dec))*floor(1/acc);
    if verbose == 1
        fprintf(['\n Found best power in Hz: power_obs = %d , and power_dec = %d ...' ...
            ' \n chi_2 prior to optimization: chi_2 = %d, \n searching phase ...'], power_obs, power_dec, chi_2_wo_phase_dec)
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval, ...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-5);
    else
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval,'TolX', 1e-5);
    end

    flag = "Rweighted2"; % select a inicial condition 'weighted','random','zeros','ones'
    Coeff_ini_obs = Coeff_Guess(0,s_p_dec, flag); %number of qbits Nq = size(flag,1) %this step may be better with separated s_p's
    Coeff_ini_dec = Coeff_Guess(0,s_p_obs, flag);
    Coeff_ini = [Coeff_ini_obs ;Coeff_ini_dec];

    if isPar== 2
        % find_phase_FSModulation_loss2q_Het(U_target,Coeff,AmpShape_obs, AmpShape_dec,s_p_obs, s_p_dec,np, ...
        %                             power_obs, power_dec, pulse_width, rho_ini, L_Spec, L_Mol, isPar) 
        % some parallelization configurations
        ms = MultiStart('FunctionTolerance',2e-15,'MaxTime', 3600,'UseParallel',true);
        problem = createOptimProblem('fmincon','x0',Coeff_ini,... 
                'objective', @(Coeff) find_phase_FSModulation_loss2q_Het(U_target,Coeff,AmpShape_obs, AmpShape_dec, ...
                                    s_p_obs, s_p_dec,np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar-1), ...
               'lb',-12/5*Coeff_ini,'ub',12/5*Coeff_ini);  %8/5Coeff_ini
        %runs the ms problem with Par_runs cpus
        Coeff_opt = run(ms,problem,Par_runs);
    elseif isPar ==1
        %function find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rho_ini, L_Spec, L_Mol)  
        %[Coeff_opt,chi_2,exitflag,output] = fminsearch(@(Coeff) ...    
        %    find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar), ...
        %    Coeff_ini,options) ;% optimization with fminsearch
        [Coeff_opt,chi_2,exitflag,output] = fminsearchbnd(@(Coeff) ...
            find_phase_FSModulation_loss2q_Het(U_target,Coeff,AmpShape_obs, AmpShape_dec, ...
                                    s_p_obs, s_p_dec,np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar), ...
            Coeff_ini,-12/5*Coeff_ini,12/5*Coeff_ini,options) ;% optimization with fminsearch with a taylored function, doesn`t work very well with upper bounds...
    else %no parallelization
        [Coeff_opt,chi_2,exitflag,output] = fminsearchbnd(@(Coeff) ...
            find_phase_FSModulation_loss2q_Het(U_target,Coeff,AmpShape_obs, AmpShape_dec, ...
                                    s_p_obs, s_p_dec,np, power_obs, power_dec,pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar), ...
            Coeff_ini,-12/5*Coeff_ini,12/5*Coeff_ini,options) ;% optimization with fminsearch with a taylored function, doesn`t work very well with upper bounds...
    end
    shape_mod_obs = make_FourierProfile([], Coeff_opt(1:s_p_obs, :), np);
    shape_mod_obs(:,2) = AmpShape_obs(:,2);
    shape_mod_dec = make_FourierProfile([], Coeff_opt(s_p_obs+1:s_p_obs +s_p_dec, :), np);
    shape_mod_dec(:,2) = AmpShape_dec(:,2);

    chi_2 = find_phase_FSModulation_loss2q_Het(U_target,Coeff_opt,shape_mod_obs, shape_mod_dec,s_p_obs, s_p_dec,np, ...
                        power_obs, power_dec, pw_obs, pw_dec, rhoini, L_Spec, L_Mol,isPar);
    if verbose == 1
        fprintf('\n found best modulation with chi^2 = %d', chi_2)        
        plot_AmpAndPhase(3, shape_mod_obs,np, pw_obs)
        plot_AmpAndPhase(4, shape_mod_dec,np, pw_dec)
    end
end

function [PowMat, chi_2_wo_phase] = find_pulse_bestpowermax_PartMultiGate_Het(PulseMod, rhoini, acc, L_Spec, L_Mol)

    rho_evol{1} = rhoini;
    rho_evol_opt{1} = rhoini;
    
    for j =1:length(PulseMod.QbitOrder)
        % unpacking PulseMod
        currentQbit = PulseMod.QbitOrder{j}{1};
        U_target = PulseMod.U_target{j};
        AmpShape_obs = PulseMod.shape_obs{j};
        AmpShape_dec = PulseMod.shape_dec{j};
        pw_obs = PulseMod.width_obs{j};
        pw_dec = PulseMod.width_dec{j};
        power_bound_obs = PulseMod.power_bound_obs{j};
        power_bound_dec = PulseMod.power_bound_dec{j};
        Par = PulseMod.Par;
    
         if currentQbit =="int" % free evolution
            pw = pw_obs;
            free_evol = ones(pw,3);
            free_evol(1:pw,1:2) = free_evol(1:pw,1:2)*1e-4;
            power = 1e-6;
            [U_gen, rho_opt] = simshapedpulse2_Par(free_evol, free_evol, pw, pw,0, 0,power,power,rho_evol{j}, L_Spec, L_Mol); 
            chi_2_opt_dec = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 ));
            power_obs = power;
            power_dec = power;
            %fprintf("\nfree evolution: \n chi_2 = %d , for pw = %d \n", chi_2_opt_dec, pw)
         elseif currentQbit == "sim"
            channel = 0;
            [power_obs, chi_2_opt] = find_pulse_bestpowermax_Het(channel, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                                    power_bound_dec, rho_evol{j}, U_target, acc);
            power_bound_obs(1,2) = power_obs;
            channel = 1;
            [power_dec, chi_2_opt] = find_pulse_bestpowermax_Het(channel, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                                    power_bound_dec, rho_evol{j}, U_target, acc);
            power_bound_dec(1,2) = power_dec;
         elseif currentQbit == "sing"
             channel = PulseMod.QbitOrder{j}{2};
             if channel ==0 
            channel = 0;
            [power_obs, chi_2_opt] = find_pulse_bestpowermax_Het(channel, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                                    power_bound_dec, rho_evol{j}, U_target, acc);
            power_dec = 1*1e-6;
             elseif channel == 1;
            [power_dec, chi_2_opt] = find_pulse_bestpowermax_Het(channel, AmpShape_obs, AmpShape_dec, pw_obs, pw_dec, power_bound_obs,...
                                                    power_bound_dec, rho_evol{j}, U_target, acc);
            power_obs = 1*1e-6;
             end
        end

        rho_evol{j+1} = U_target*rho_evol{j}*U_target';

        chi_2_wo_phase(j) = chi_2_opt;

        PowMat(1,j) = power_obs;
        PowMat(2,j) = power_dec;
    end
end

function [power, chi_2_opt] = find_pulse_bestpowermax_Het(channel, shape_obs, shape_dec, pw_obs, pw_dec, power_bound_obs,...
                                            power_bound_dec, rho_ini, U_target, acc) 
    % from a pulse_width especification and a fixed regular shape modulation (e.g. hard.RF, gauss.RF, isech.RF)
    % finds the best amplitude to do a pi pulse.
    % this code doesn't run the actual maxpower, since it is assumed that it is experimentally forbiden
    % default should be for 180 degrees on first qbit
    arguments 
        channel double % is 0 for observer and 1 for decoupler
        shape_obs (:,3) % loaded from an '.RF' file
        shape_dec (:,3) % loaded from an '.RF' file
        pw_obs double
        pw_dec double
        power_bound_obs (1,2) double = [0 2*1e4] %works for 2*1eN
        power_bound_dec (1,2) double = [0 2*1e4] %works for 2*1eN
        rho_ini (:,:) double = eye(8)/8 % should be density matrix, usually in state |0^N>
        U_target (:,:) double = kron([0 1; 1 0], kron(eye(2), eye(2)))
        acc double = 1e-4 %must be powers of 10
    end
    if channel ==0
        minpower = power_bound_obs(1,1);
        maxpower = power_bound_obs(1,2);
        power_dec = power_bound_dec(1,2);
    elseif channel ==1
        minpower = power_bound_dec(1,1);
        maxpower = power_bound_dec(1,2);
        power_obs = power_bound_obs(1,2);
    end

    % routine for non rounded maxpower*1e-N
    lexp = floor(log10(maxpower));
    hexp = floor(maxpower*10^(-lexp+1));
    hexp = hexp/10;
      
    %setting initial values
    Npulse = floor(1/acc);
    power = maxpower;
    chi_2_opt = 100;
    
    % 
    %U_target = kron(sx,kron(id,id));

    if Npulse>10
        [power, chi_2_opt] = find_pulse_bestpowermax_Het(channel, shape_obs, shape_dec, pw_obs, pw_dec, ...
                                            power_bound_obs, power_bound_dec,rho_ini, U_target,acc*10); 
        Nlim = 2*10 -1;
        maxpower_loop = power;
        if (power)<=maxpower*(1-acc/2) %some conditions to not exceed maximum power bound.
            maxpower_loop = maxpower_loop +maxpower*acc*10/2;
        end
        for r=1:Nlim % runs +-10% maxpower over last found power
            %[U, roh] = simshapedpulse2_Par(shape_obs,shape_dec,pw_obs,pw_dec,... 
            %                    ramp_obs,ramp_dec,power_obs,power_dec,roha
            if channel == 0
                [U_gen, ~] = simshapedpulse2_noPar(shape_obs,shape_dec,pw_obs,pw_dec,... 
                                    0,0,(maxpower_loop -r*acc*1e4),power_dec,rho_ini); 
            elseif channel ==1
                [U_gen, ~] = simshapedpulse2_noPar(shape_obs,shape_dec,pw_obs,pw_dec,... 
                                        0,0,power_obs,(maxpower_loop -r*acc*1e4),rho_ini); 
            end
            chi_2_r = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 ));
            if chi_2_r<chi_2_opt
                p = maxpower_loop -r*acc*1e4;
                if p >= minpower
                    power = p;
                    chi_2_opt = chi_2_r;
                    if chi_2_opt < acc
                        return
                    end
                end
            end
        end
    else
        for r=1:hexp*Npulse 
            r=r/2;
            %[U, roh] = simshapedpulse2_Par(shape_obs,shape_dec,pw_obs,pw_dec,... 
            %                    ramp_obs,ramp_dec,power_obs,power_dec,roha
            if channel == 0
                [U_gen, ~] = simshapedpulse2_noPar(shape_obs,shape_dec,pw_obs,pw_dec,... 
                                    0,0,(1-r*acc)*maxpower,power_dec,rho_ini); 
            elseif channel ==1
                [U_gen, ~] = simshapedpulse2_noPar(shape_obs,shape_dec,pw_obs,pw_dec,... 
                                        0,0,power_obs,(1-r*acc)*maxpower,rho_ini); 
            end
            chi_2_r = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 ));
            if chi_2_r<chi_2_opt
                p = (1-r*acc)*maxpower;
                if p >= minpower
                    power = p;
                    chi_2_opt = chi_2_r;
                    if chi_2_opt < acc
                        return
                    end
                end
            end
        end
    end
end
   
function [PulseMod, chi_2_fullFS, rho_evol] = get_FSMod_PartMultigate_Het(PulseMod,rhoini, L_Spec,L_Mol, verbose)
    % this function partitions the optimization by the Gate order given in structured form
    rho_evol{1} = rhoini;
    rho_evol_opt{1} = rhoini;
    acc = 1e-4;

    for j =1:length(PulseMod.QbitOrder)
        if verbose ==1
            fprintf('\noptimizing gate %d of %d', j , length(PulseMod.QbitOrder))
        end
        % unpacking PulseMod
        currentQbit = PulseMod.QbitOrder{j};
        U_target = PulseMod.U_target{j};
        pw_obs = PulseMod.width_obs{j};
        pw_dec = PulseMod.width_dec{j};
        power_bound_obs = PulseMod.power_bound_obs{j};
        power_bound_dec = PulseMod.power_bound_dec{j};
        power = power_bound_obs(2);
        s_pVec(1) = PulseMod.s_p_obs{j};
        s_pVec(2) = PulseMod.s_p_dec{j};
        s_aVec(1) = PulseMod.s_a_obs{j};
        s_aVec(2) = PulseMod.s_a_dec{j};
        np = PulseMod.np_obs{j};
        Par = PulseMod.Par;

        switch currentQbit{1}
            case {"int", 'grad'} % free evolution or gradient
                pw = pw_obs;
                free_evol = ones(pw,3);
                free_evol(1:pw,1:2) = free_evol(1:pw,1:2)*1e-4;
                power = 1e-6;
                if currentQbit{1} =="int" %free evolution
                    [U_gen, rho_opt] = simshapedpulse2_Par(free_evol, free_evol, pw, pw,0, 0,power,power,rho_evol{j}, L_Spec, L_Mol); 
                    chi_2 = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 )) + abs(angle(U_gen(1,1)) - angle(U_target(1,1)));
                    if verbose ==1
                        fprintf("\nModulaton metrics (free evolution): \n chi_2 = %d , for pw = %d \n", chi_2, pw)
                    end
                    rho_evol{j+1} = U_target*rho_evol{j}*U_target';
                else % gradient                    
                    if verbose == 1
                        fprintf('running gradient on gate number %d', j)
                    end
                    D = diag(rho_evol{j});
                    rho1 = eye(4);
                    for k2 =1:4
                        rho1(k2,k2) = D(k2);
                    end
                    rho_evol{j+1} = rho1;
                    chi_2 =0;
                end
            shape_mod_obs = free_evol;
            shape_mod_dec = free_evol;
            % actual gate
            case{'sim'}
                [shape_mod_obs, shape_mod_dec,chi_2] = ...
                    get_FSMod_2q_Het(s_aVec, s_pVec, np, pw_obs, pw_dec, power_bound_obs,power_bound_dec, U_target, rho_evol{j}, L_Spec,L_Mol, verbose, Par);
                rho_evol{j+1} = U_target*rho_evol{j}*U_target';
            case{'sing'}
                chn = currentQbit{2};
                if chn == 0
                 [shape_mod_obs,chi_2] = get_FSMod_1q_Het(s_aVec(1), s_pVec(1),...
                        np, pw_obs,power, U_target, rho_evol{1}, L_Spec,L_Mol, verbose, Par, 0);
                 shape_mod_dec = shape_mod_obs*0;
                elseif chn == 1
                 [shape_mod_dec,chi_2] = get_FSMod_1q_Het(s_aVec(2), s_pVec(2),...
                        np, pw_dec,power, U_target, rho_evol{1}, L_Spec,L_Mol, verbose, Par, 1);
                 shape_mod_obs = shape_mod_dec*0;
                end
                rho_evol{j+1} = U_target*rho_evol{j}*U_target';
        end        
             
        
        chi_2_fullFS(j) = chi_2;
       
        PulseMod.shape_obs{j} = shape_mod_obs;
        PulseMod.shape_dec{j} = shape_mod_dec;
        
    end

end

function [shape_mod_obs, shape_mod_dec, chi_2] = get_FSMod_2q_Het(s_aVec, s_pVec, np, ... 
                     pulse_width_obs, pulse_width_dec, power_bound_obs,power_bound_dec, U_target, rhoini, L_Spec,L_Mol, verbose, Par)
    %paralellization runs as follows: isPar = 0 for non-parallel, 1 for parallelization defined outside the function as a parfor
    % 2 for parallelization inside as a multistart, it then has a second argument as Parruns = number of starting points
    arguments
        s_aVec (1,:) % vector with number of coefficients in amplitude sum
        s_pVec (1,:) % vector with number of coefficients in phase sum [s_p_Qbit1 s_p_Qbit2]
        np double % number of sample points in time
        pulse_width_obs double % time of the pulse
        pulse_width_dec double % time of the pulse
        power_bound_obs (1,:) % power in obs channel
        power_bound_dec (1,:) % power in both channels
        U_target (:,:) % Target unitary
        rhoini (:,:) % inicial state
        L_Spec (:,:,:) % Spectropar variable 
        L_Mol (:,:,:)
        verbose double
        Par (1,2)
    end
    % unpacking variables
    isPar = Par(1);
    Par_runs = Par(2);
    s_p_obs = s_pVec(1);
    s_p_dec = s_pVec(2);
    s_a_obs = s_aVec(1);
    s_a_dec = s_aVec(2);
    sum_sa = s_a_obs +s_a_dec;
    sum_sp = s_p_obs+ s_p_dec; 
    power_obs = power_bound_obs(2);
    power_dec = power_bound_dec(2);
    % for now the pulses have the same width and maximum power
    pulse_width = pulse_width_obs; 

    % setting optimization options
    Maxeval = (sum_sa+sum_sp)*5*1e2; %*20/5
    if verbose ==1 
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval, ...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-5);
        if isPar == 1
            fprintf(['\n searching for best fourier coefficients with \n observer Sa Sp = %d , %d ...' ...
                '\n decoupler Sa, Sp = %d, %d ... \n and with Parallel runs determined with parfor...'],...
                s_a_obs, s_p_obs, s_a_dec, s_p_dec)
        elseif isPar == 2
            fprintf(['\n searching for best fourier coefficients with \n observer Sa Sp = %d , %d ...' ...
                '\n decoupler Sa, Sp = %d, %d ... \n with Parallel %d runs ...'], s_a_obs, s_p_obs, s_a_dec, s_p_dec, Par_runs)
        end
    else
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval,'TolX', 1e-5);
    end

    % % Modulating 2 qubit in a heteronuclear system
    flag = ["Rweighted2"];%,"Rweighted2"] ;
    Coeff_ini_obs = Coeff_Guess(s_a_obs,s_p_obs, flag); %number of qbits Nq = size(flag,1) %this step may be better with separated s_p's
    Coeff_ini_dec = Coeff_Guess(s_a_dec,s_p_dec, flag);
    Coeff_ini = [Coeff_ini_obs ;Coeff_ini_dec];
    %Coeff_ini = Coeff_Guess(s_a,s_p, flag); %number of qbits Nq = size(flag,1)
    Coeff_ini_LB = 0*Coeff_ini;
    Coeff_ini_LB(s_a_obs+1:s_a_obs + s_p_obs, :) = -1*Coeff_ini(s_a_obs+1:s_a_obs + s_p_obs, :);
    Coeff_ini_LB(sum_sa+s_p_obs+1:sum_sa+sum_sp,:)= -1*Coeff_ini(sum_sa+s_p_obs+1:sum_sa+sum_sp,:);

    if isPar== 2
        ms = MultiStart('FunctionTolerance',2e-15,'MaxTime', 3600,'UseParallel',true);
        % for 2 frequencies (2 qbits)
           problem = createOptimProblem('fmincon','x0',Coeff_ini,...
            'objective', @(Coeff) find_FSModulation_loss2q(U_target,Coeff,s_aVec,s_pVec,np, ...
            power_obs,power_dec, pulse_width, rhoini, L_Spec, L_Mol, (isPar-1)), 'lb',Coeff_ini_LB,'ub',12/5*Coeff_ini); 
        %runs the ms problem with Par_runs cpus
        Coeff_opt = run(ms,problem,Par_runs);
    elseif isPar ==1
        % optimization with fminsearch with parrallelization defined outside the function
        [Coeff_opt,chi_2,~,output] = fminsearch(@(Coeff) ...
            find_FSModulation_loss2q(U_target,Coeff,s_aVec,s_pVec,np, power_obs,power_dec, ...
            pulse_width, rhoini, L_Spec, L_Mol, isPar), Coeff_ini, options) ;
    elseif isPar ==0
        % optimization with fminsearch without parallelization
        [Coeff_opt,chi_2,~,output] = fminsearch(@(Coeff) ...
            find_FSModulation_loss2q(U_target,Coeff,s_aVec,s_pVec,np, power_obs,power_dec, ...
            pulse_width, rhoini, L_Spec, L_Mol, isPar),Coeff_ini, options) ;
    end
    
    
    % Amp and Phase Coefficients spliting 
    A_Coef_obs = Coeff_opt(1:s_a_obs, :);
    P_Coef_obs = Coeff_opt(s_a_obs+1:s_a_obs + s_p_obs, :);
    [p1,~] = size(Coeff_opt);
    if (p1 > s_a_obs+s_p_obs) %split obs/dec
        A_Coef_dec = Coeff_opt(s_a_obs + s_p_obs + 1 : sum_sa + s_p_obs, :);
        P_Coef_dec = Coeff_opt(sum_sa + s_p_obs + 1 : sum_sa +sum_sp, :);
        shape_mod_dec = make_FourierProfile(A_Coef_dec, P_Coef_dec, np);
    end
    shape_mod_obs = make_FourierProfile(A_Coef_obs, P_Coef_obs, np);
    
    % getting the final [U, rhofin]
    %function [U, roh] = simshapedpulse2_noPar(shape_obs,shape_dec,pw_obs,pw_dec,... 
    %ramp_obs,ramp_dec,power_obs,power_dec,roha, L_Spec, L_Mol)
    [U_opt, ~] = simshapedpulse2_noPar(shape_mod_obs,shape_mod_dec,pulse_width,pulse_width, ...
                    0,0,power_obs,power_dec,rhoini);% , L_Spec, L_Mol);
    
    % Testing metrics on Outputs
    %abs(U_opt - U_target);
    [f1,~] = size(U_target);
    F_jhon = 1-abs(trace(U_opt'*U_target))/f1;
    chi_2 = sum(sum( real(U_opt - U_target).^2 + imag(U_opt - U_target).^2 ));
    
    if verbose ==1
        % Plotting Routines, to see pulse shapes
        plot_AmpAndPhase(1,shape_mod_obs, np, pulse_width);
        plot_AmpAndPhase(2,shape_mod_dec, np, pulse_width);
        fprintf('\n Modulation metrics: \n chi_2 = %d \n Fid_John = %d', chi_2, F_jhon)
    end
end

function [shape_mod_obs, chi_2] = get_FSMod_1q_Het(s_a, s_p, np, ... 
                     pulse_width, power, U_target, rhoini, L_Spec,L_Mol, verbose, Par, channel)
    %paralellization runs as follows: isPar = 0 for non-parallel, 1 for parallelization defined outside the function as a parfor
    % 2 for parallelization inside as a multistart, it then has a second argument as Parruns = number of starting points
    arguments
        s_a double % vector with number of coefficients in amplitude sum
        s_p double % vector with number of coefficients in phase sum [s_p_Qbit1 s_p_Qbit2]
        np double % number of sample points in time
        pulse_width double % time of the pulse
        power double % power in both channels
        U_target (:,:) % Target unitary
        rhoini (:,:) % inicial state
        L_Spec (:,:,:) % Spectropar variable 
        L_Mol (:,:,:)
        verbose double
        Par (1,2)
        channel double % 1 for obs , 2 for dec
    end
    isPar = Par(1);
    Par_runs = Par(2);
    
    % for now the pulses have the same width and maximum power
    
    % setting optimization options
    Maxeval = (s_a+s_p)*4*1e3;
    
    if verbose ==1 
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval, ...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-5);
        if isPar == 1
            fprintf(['\n searching for best fourier coefficients with \n observer Sa Sp = %d , %d ...' ...
                '\n and with Parallel runs determined with parfor...'],...
                s_a, s_p)
        elseif isPar == 2
            fprintf(['\n searching for best fourier coefficients with \n observer Sa Sp = %d , %d ...' ...
                '\n with Parallel %d runs ...'], s_a, s_p, Par_runs)
        end
    else
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval,'TolX', 1e-5);
    end
    

    % % Modulating 2 qubit in a heteronuclear system
    flag = ["Rweighted2"];%,"Rweighted2"] ;
    Coeff_ini= Coeff_Guess(s_a,s_p, flag); %number of qbits Nq = size(flag,1) %this step may be better with separated s_p's
    %Coeff_ini = Coeff_Guess(s_a,s_p, flag); %number of qbits Nq = size(flag,1)
    Coeff_ini_LB = 0*Coeff_ini;
    Coeff_ini_LB(s_a+1:s_a+ s_p, :) = -1*Coeff_ini(s_a+1:s_a+ s_p, :);
    %Coeff_ini_LB(s_a+s_p+1:sum_sa+sum_sp,:)= -1*Coeff_ini(sum_sa+s_p_obs+1:sum_sa+sum_sp,:);

    if isPar== 2
        ms = MultiStart('FunctionTolerance',2e-15,'MaxTime', 3600,'UseParallel',true);
        % for 2 frequencies (2 qbits)
        %find_FSModulation_loss1q(U_target,Coeff,s_a,s_p,np, power, pulse_width, rho_ini) 
           problem = createOptimProblem('fmincon','x0',Coeff_ini,...
            'objective', @(Coeff) find_FSModulation_loss1q_Het(U_target,Coeff,s_a,s_p,np, ...
            power, pulse_width, rhoini, L_Spec, L_Mol, (isPar-1),channel), 'lb',Coeff_ini_LB,'ub',12/5*Coeff_ini); 
        %runs the ms problem with Par_runs cpus
        Coeff_opt = run(ms,problem,Par_runs);
    elseif isPar ==1
        % optimization with fminsearch with parrallelization defined outside the function
        [Coeff_opt,chi_2,~,output] = fminsearch(@(Coeff) ...
            find_FSModulation_loss1q_Het(U_target,Coeff,s_a,s_p,np, power, ...
            pulse_width, rhoini, L_Spec, L_Mol, isPar,channel), Coeff_ini, options) ;
    elseif isPar ==0
        % optimization with fminsearch without parallelization
        [Coeff_opt,chi_2,~,output] = fminsearch(@(Coeff) ...
            find_FSModulation_loss1q_Het(U_target,Coeff,s_a,s_p,np,power, ...
            pulse_width, rhoini, L_Spec, L_Mol, isPar,channel),Coeff_ini, options) ;
    end
    %}
    
    % Amp and Phase Coefficients spliting 
    A_Coef_obs = Coeff_opt(1:s_a, :);
    P_Coef_obs = Coeff_opt(s_a+1:s_a+ s_p, :);
    [p1,~] = size(Coeff_opt);
    shape_mod_obs = make_FourierProfile(A_Coef_obs, P_Coef_obs, np);
    
    % getting the final [U, rhofin]
    %function [U, roh] = simshapedpulse2_noPar(shape_obs,shape_dec,pw_obs,pw_dec,... 
    %ramp_obs,ramp_dec,power_obs,power_dec,roha, L_Spec, L_Mol)
    [U_opt, ~] = shapedpulse_Par(shape_mod_obs,pulse_width, [],0,power,rhoini, L_Spec, L_Mol);
    
    % Testing metrics on Outputs
    %abs(U_opt - U_target);
    [f1,~] = size(U_target);
    F_jhon = 1-abs(trace(U_opt'*U_target))/f1;
    chi_2 = sum(sum( real(U_opt - U_target).^2 + imag(U_opt - U_target).^2 ));
    
    if verbose ==1
        % Plotting Routines, to see pulse shapes
        plot_AmpAndPhase(1,shape_mod_obs, np, pulse_width);
        %plot_AmpAndPhase(2,shape_mod_dec, np, pulse_width);
        fprintf('\n Modulation metrics: \n chi_2 = %d \n Fid_John = %d', chi_2, F_jhon)
    end
end

function [PulseMod , chi_2, rho_evol_opt] = get_FSMod_PartMultigate_Hom(PulseMod, rhoini, L_Spec, verbose);
%finds the best modulations for amplitude and phaseon each unitary applied on a circuit, assumes that the Pulse
% shapes and properties are all in PulseMod structure and that unitaries are partitioned 
% and ordered in PulseMod.U_target. 
arguments
    PulseMod struct %structure with pulse properties
    rhoini (:,:) %inicial state
    L_Spec (1,:) %spectropar as local variable %L_Mol (:,:,:,:) % Molecules Ix's and Iy's stored as a local variable
    verbose double % 0 for False, 1 for true
end
    rho_evol{1} = rhoini;
    rho_evol_opt{1} = rhoini;
    acc = 1e-4;

    for j =1:length(PulseMod.QbitOrder)
        if verbose ==1
            fprintf('\noptimizing gate %d of %d', j , length(PulseMod.QbitOrder))
        end
        % unpacking PulseMod
        currentQbit = PulseMod.QbitOrder{j};
        U_target = PulseMod.U_target{j};
        pw_obs = PulseMod.width_obs{j};
        power_bound_obs = PulseMod.power_bound_obs{j};
        s_a = PulseMod.s_a_obs{j};
        s_p = PulseMod.s_p_obs{j};
        np = PulseMod.np_obs{j};
        Par = PulseMod.Par;

        switch currentQbit{1}
            case {"int", 'grad'} % free evolution or gradient
                pw = pw_obs;
                free_evol = ones(pw,3);
                free_evol(1:pw,1:2) = free_evol(1:pw,1:2)*1e-4;
                power = 1e-6;
                if currentQbit{1} == "int" %free evolution
                    %shapedpulse_Hom_Par(shape_mod,pulse_width,shape_mod(:,1),0,power,rho_ini, L_Spec, L_Mol); 
                        [U_gen, rho_opt] = shapedpulse_Hom_Par(free_evol, pw, [],0,power,rho_evol{j}, L_Spec, L_Mol); 
                        chi_2 = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 )) + abs(angle(U_gen(1,1)) - angle(U_target(1,1)));
                        if verbose ==1
                            fprintf("\nModulaton metrics (free evolution): \n chi_2 = %d , for pw = %d \n", chi_2, pw)
                        end
                        rho_evol{j+1} = U_target*rho_evol{j}*U_target';
                else % gradient                    
                    if verbose == 1
                        fprintf('running gradient on gate number %d', j)
                    end
                    D = diag(rho_evol{j});
                    rho1 = eye(4);
                    for k2 =1:4
                        rho1(k2,k2) = D(k2);
                    end
                    rho_evol{j+1} = rho1;
                    chi_2 =0;
                end
                shape_mod_obs = free_evol;
            case{'sing'}             

                % setting resonance frequency
                target_qbit = currentQbit{2};
                offset = PulseMod.offset{j};
                L_Mol = focusfrequency(target_qbit, offset);      

                %running pulse optimization
                [shape_mod_obs,chi_2(j)] = get_FSMod_1q_Hom(s_a, s_p, np, pw_obs, ...
                    power_bound_obs(2), U_target, rho_evol{j}, L_Spec,L_Mol, verbose, Par);

                rho_evol{j+1} = U_target*rho_evol{j}*U_target';
        end
        % finds full modulation for amplitude and phase
        %[shape_mod_obs,chi_2(j)] = get_FSMod_1q_Hom(s_a, s_p, np, pulse_width_obs, ...
         %   power_bound_obs(2), U_target, rho_evol{j}, L_Spec,L_Mol, verbose, Par);
          
        % getting the final [U, rhofin]
        %[U_gen, rhofin] = shapedpulse_Hom_Par(shape,pw,[],0,(maxpower_loop -r*acc*1e4),rho_ini, L_Spec, L_Mol); 
        [U_opt, ~] = shapedpulse_Hom_Par(shape_mod_obs,pw_obs,[],0,power_bound_obs(2),rho_evol{j}, L_Spec, L_Mol);
        %chi_22(j) = sum(sum( real(U_opt - U_target).^2 + imag(U_opt - U_target).^2 ));
        rho_evol{j+1} = U_target*rho_evol{j}*U_target';
        rho_evol_opt{j+1} = U_opt*rho_evol_opt{j}*U_opt';
        PulseMod.shape_obs{j} = shape_mod_obs;
    end


end

function L_Mol = focusfrequency(target_qbit, offset)
   % rewrites the L_Mol but with the desired resonance frequency of target qbit + specified offset;
    global mol spectro
        str = append(mol.nome, num2str(target_qbit+1));
        funcHandle = str2func(str);
        funcHandle();
        mol.dq = mol.dq - offset;
        [L_Spectro, L_Mol] = fetch_global_val(spectro,mol);
        [L_Mol(:,:,2) L_Mol(:,:,1)] = GenH(mol);
end


function L_Mol = focusfrequency_Het(target_qbit, offset)
   % rewrites the L_Mol but with the desired resonance frequency of target qbit + specified offset;
    global mol spectro
        str = mol.nome;
        funcHandle = str2func(str);
        funcHandle();
        mol.dq(target_qbit) = mol.dq(target_qbit) - offset;
        [L_Spectro, L_Mol] = fetch_global_val(spectro,mol);
        [L_Mol(:,:,2) L_Mol(:,:,1)] = GenH(mol);
end

function [shape_mod_obs,chi_2] = get_FSMod_1q_Hom(s_a, s_p, np, pulse_width, power, U_target, rhoini, L_Spec,L_Mol, verbose, Par)
    %given the parameters of a experimental gate, and number of coefficients for a Fourier Pulse, returns the shape for both
    %amplitude and phase as a npx3 [phase amp time_steps] matrix. 
    %obs: L_Spec and L_Mol are only used when parallelization is on
    arguments
        s_a double = 3 % number of coeff to sum in amplitude
        s_p double  = 5 % number of coeff to sum in phase
        np double = 200 % number of points
        pulse_width double = 200 % time duration of pulse in (us)
        power double = 20*1e3 % max_power on channel
        U_target (:,:) = kron([0 1; 1 0], eye(2)) %unitary evolution target
        rhoini (:,:) = [1 0; 0 0] % inicial state
        L_Spec (1,:) =  [10.0000   24.8750] 
        L_Mol (:,:,:) = zeros(8,8,6)
        verbose double = 0 %1 for verbose
        Par (1,:) = [0 44] % Par = [isPar Par_runs] decides if Parallelization is used, 2 for par in function, 1 outside, 0 unpar
    end
    isPar = Par(1);
    Par_runs = Par(2);

    % setting optimization options
    Maxeval = (s_a+s_p)*4*1e3;
    if verbose ==1 
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval, ...
            'PlotFcns', {@optimplotx, @optimplotfunccount, @optimplotfval} ,'TolX', 1e-5);
        if isPar == 1
            fprintf('\n searching for best fourier coefficients with s_a s_p = %d , %d and with Parallel runs determined with parfor...', s_a, s_p)
        elseif isPar == 2
            fprintf('\n searching for best fourier coefficients with s_a s_p = %d , %d and with %d Parallel runs...', s_a, s_p, Par_runs)
        end
    else
        options = optimset('MaxFunEvals', Maxeval,'MaxIter', Maxeval,'TolX', 1e-5);
    end

    % % Modulating 1 qubit in a homonuclear system
    flag = ["Rweighted2"] ;
    Coeff_ini = Coeff_Guess(s_a,s_p, flag); %number of qbits Nq = size(flag,1)
    Coeff_ini_LB = 0*Coeff_ini;
    Coeff_ini_LB(s_a+1:s_a + s_p, :) = -1*Coeff_ini(s_a+1:s_a + s_p, :);
    
    if isPar== 2
        [p1,q1] = size(Coeff_ini);
        % some parallelization configurations
        ms = MultiStart('FunctionTolerance',2e-15,'MaxTime', 3600,'UseParallel',true);
        problem = createOptimProblem('fmincon','x0',Coeff_ini,... 
                'objective', @(Coeff) find_FSModulation_loss1q_Hom(U_target,Coeff,s_a,s_p,np, power,pulse_width, rhoini, L_Spec, L_Mol,(isPar-1)), ...
               'lb',12/5*Coeff_ini_LB,'ub',12/5*Coeff_ini);  %8/5Coeff_ini
        %runs the ms problem with Par_runs cpus
        Coeff_opt = run(ms,problem,Par_runs);
    elseif isPar ==1
        [Coeff_opt,chi_2,exitflag,output] = fminsearch(@(Coeff) ...
         find_FSModulation_loss1q_Hom(U_target,Coeff,s_a,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar), Coeff_ini, options) ;
            % optimization with fminsearch still parallel but parallelization is outside the function
    elseif isPar ==0
        [Coeff_opt,chi_2,exitflag,output] = fminsearch(@(Coeff) ...
         find_FSModulation_loss1q_Hom(U_target,Coeff,s_a,s_p,np, power, pulse_width, rhoini, L_Spec, L_Mol,isPar), Coeff_ini, options) ;
            % optimization with fminsearch not parallel
    end

     % Amp and Phase Coefficients spliting 
    A_Coef_obs = Coeff_opt(1:s_a, :);
    P_Coef_obs = Coeff_opt(s_a+1:s_a + s_p, :);
    shape_mod_dec = [0,0,0]; % uncomment for single channel pulses
    [p1,~] = size(Coeff_opt);
    if (p1 > s_a+s_p) %split obs/dec
        A_Coef_dec = Coeff_opt(s_a + s_p + 1 : 2*s_a + s_p, :);
        P_Coef_dec = Coeff_opt(2*s_a + s_p + 1 : 2*(s_a + s_p), :);
        shape_mod_dec = make_FourierProfile(A_Coef_dec, P_Coef_dec, np);
    end
    shape_mod_obs = make_FourierProfile(A_Coef_obs, P_Coef_obs, np);
       
    % getting the final [U, rhofin]
    %[U_opt, rohfin] = shapedpulse_Hom(shape_mod_obs,pulse_width,shape_mod_obs(:,1),ramp_obs, power,rhoini) ; %,L_Spec, L_Mol);
    [U_opt, rhofin] = shapedpulse_Hom_Par(shape_mod_obs,pulse_width,shape_mod_obs(:,1),0,power,rhoini, L_Spec, L_Mol);
        
    % Testing metrics on Outputs
    %abs(U_opt - U_target);
    [f1,~] = size(U_target);
    F_jhon = 1-abs(trace(U_opt'*U_target))/f1;
    %chi_2 = sum(sum(abs(abs(U_opt) - abs(U_target))));
    chi_2 = sum(sum( real(U_opt - U_target).^2 + imag(U_opt - U_target).^2 ));
    
    if verbose ==1
        % Plotting Routines, to see pulse shapes
        plot_AmpAndPhase(1,shape_mod_obs, np, pulse_width);
        %plot_AmpAndPhase(2,shape_mod_dec, np, pulse_width);
        fprintf('\n Modulation metrics: \n chi_2 = %d \n Fid_John = %d', chi_2, F_jhon)
    end

end


function result = find_modulation_loss1q(U_target,chute,rho_ini,flag) 
    % Given a gaussian pulse modulation, this function finds the best time to stop
    % it, so it equals a U_targert.
    % chute = inicial time guess
    % U_taget unitary (on time) target is 4x4 has both qbits
 
    shape_gauss = make_gassian_profile(250);
    
    % >>> Modulated pulse evolution <<<<<
    % pulse shape is a matrix = (phase, amplitude, timesteps)
    %[U roh] = shapedpulse2(shape,pw,phase,ramp,power,roha);
    [U, rhofin] = shapedpulse2(shape_gauss,chute,90,0,20000,rho_ini);
    
    % Function to be minimized chi^2 
    chi_2 = sum(sum( real(U - U_target).^2 + imag(U - U_target).^2 ));
    
    % Two outputs: if flag is given return U, otherwise return chi_2   
    if nargin > 3 % nargin -> numero de argumentos da funcao 
        result = U; 
    else
        result = chi_2;
    end

end

function [U_evol, rho_evol] = rho_in_timeseries(shape_obs, shape_dec, pw_obs, pw_dec, ramp_obs, ramp_dec, ...
    power_obs, power_dec,  channel_num, rhoini, np , pulsetype)
    % this function gets a pulse and a inicial state and gives the time evolution of the state until the final time pw
    % it works for multiple types of gates, either 1 or 2 channels, and bot homonuclear and heteronuclear systems
    
    %------------- has to have Parallelization options later --------%

    % DEFAULT: run obs and dec are the same, even when applying single channel pulses.

    % inputs shape_obs or shape_dec are either saved .RF, i.e. 'char' files with (phase, amplitude, timestep) or matrices
    % shape_obs/_dec = either  '.RF' file or a 3xN matrix with (amplitude, phase, timestep)
    % pw_obs/dec = float : pulse width; ramp_obs/dec = 1xnp vector : ramp on phase ; power_obs/dec = float :  max power of channel
    % channel_num = int : number of channels, usually 1 or 2. ; rhoini = matrix: inicial state, density matrix ; np = integer: number of points (has to be
    % les than pw_obs/dec; pulsetype = string : can be either 'shapedpulse2','simshapedpulse2_noPar', 'shapedpulse_Hom' , 
    % OBS: '_Hom' suffix should be used for homonuclear pulses
    rho_evol{1} = rhoini;
    [p1,~] = size(rhoini);
    U_evol{1} = eye(p1);

    %changes order of columns if called a saved RF file <--- this should change when order is fixed
    if isa(shape_obs,"char")
        %shape_obs = RF_file2Matrix(shape_obs);
        shape_obs = load(shape_obs);
        %phase_obs = shape_obs(2,:);
    end
    if isa(shape_dec,'char')
        shape_dec = load(shape_dec);
        %phase_dec = shape_dec(2,:);
    end
    
    % defining pulse
    if channel_num == 2
        pulse_func = @simshapedpulse2_noPar; %should later be adapted to have simshapedpulse_Hom
    elseif channel_num == 1 
        if strcmp(pulsetype,'shapedpulse2')
            pulse_func = @shapedpulse2;
        elseif strcmp(pulsetype,'shapedpulse_Hom')
            pulse_func = @shapedpulse_Hom;
        end
    end
    
    % setting timesteps
    dtau_obs = pw_obs/np;
    dtau_dec = pw_dec/np;
    tau_obs =dtau_obs;
    tau_dec = dtau_dec;

    % loop  for saving each dtau in rho_evol
    U = 1; 
    for j = 1:np
        if channel_num == 2
            [U, rhofin] = pulse_func(shape_obs,shape_dec,tau_obs,tau_dec, ramp_obs,ramp_dec,power_obs,power_dec,rhoini);
            atau_dec(j) = tau_dec;
            tau_dec = tau_dec + dtau_dec;
        elseif channel_num == 1
            [U, rhofin] = pulse_func(shape_obs,tau_obs,shape_obs(:,1),ramp_obs,power_obs,rhoini);
        end
        rho_evol{j+1} = rhofin;
        U_evol{j+1} = U;
        atau_obs(j) = tau_obs;
        tau_obs = tau_obs + dtau_obs;
    end
   clearvars pulse_func 

end

function [U_opt, rho_opt] = get_UAndRho_Het(PulseMod, rhoini,PowMatrix, L_Spec,L_Mol, isfullmod)
 %once optimization is done, run this function to get the U_opt and rho_opt
    rho_opt{1} = rhoini;
    for j =1:length(PulseMod.QbitOrder)
        currentQbit = PulseMod.QbitOrder{j};
        AmpShape_obs = PulseMod.shape_obs{j};
        AmpShape_dec = PulseMod.shape_dec{j};
        pw_obs = PulseMod.width_obs{j};
        pw_dec = PulseMod.width_dec{j};
        if (isfullmod)
            power_bound_obs = PulseMod.power_bound_obs{j};
            power_obs = power_bound_obs(2);
            power_bound_dec = PulseMod.power_bound_dec{j};
            power_dec = power_bound_dec(2);
        else
            power_obs = PowMatrix(1,j);
            power_dec = PowMatrix(2,j);
        end
        switch currentQbit
            case {0, 'grad'} % free evolution or gradient
                pw = pw_obs;
                free_evol = ones(pw,3);
                free_evol(1:pw,1:2) = free_evol(1:pw,1:2)*1e-4;
                power = 1e-6;
                if currentQbit ==0 %free evolution
                    [U_opt{j}, rho_opt{j+1}] = simshapedpulse2_Par(free_evol, free_evol, pw, pw,0, 0,power,power,rho_opt{j}, L_Spec, L_Mol); 
                else % gradient                    
                    D = diag(rho_opt{j});
                    rho1 = eye(4);
                    for k2 =1:4
                        rho1(k2,k2) = D(k2);
                    end
                    rho_opt{j+1} = rho1;
                end
            case {1} % actual gate
            [U_opt{j}, rho_opt{j+1}] = simshapedpulse2_Par(AmpShape_obs, AmpShape_dec,pw_obs, pw_dec, 0, 0, ...
                            power_obs,power_dec,rho_opt{j}, L_Spec, L_Mol); 
            rho_opt{j+1} = U_opt{j}*rho_opt{j}*U_opt{j}';
         end        
     end
end


%{
function shape = RF_file2Matrix(RF_file)
    %gets RF file and changes it to a matrix, also changes the order of the two first columns,%%%  <---- this should be changed later
    shp_2 = load(RF_file);
    shp = zeros(length(shp_2),3);
    shp(:,1) = shp_2(:,2);
    shp(:,2) = shp_2(:,1);
    shp(:,3) = ones(length(shp_2),1);
    shape = shp;
end
%}

function pulse_shape = make_gassian_profile(np)
    % returns a pulse as a gaussian in np points
    % Reproduces the Gaussian Pulse adapted to Varian
    
    % >>>>> Defining the pulse profile as a matrix [phase amplitude time_partition] <<<<< 
    % Gaussian Pulse:  250 points to reproduce the Varian archives
    % This pulse is amplitude modulated and selectively excites 
    % a bandwidth (Hz) approximately equal to 2e+6/pulse_length (us).
    % gauss: Amplitude modulated gaussian pulse: b1(x)=1023*exp(-x^2)

    %np = 250; default value
    x = -np/2;
    dx = 1;
    
    phase = 0.0; % constant phase
    time_part = 1.0;
    
    %std_dev = pow2(ceil(log2(np/2))); % finds highest power of 2 that aproximates np/2
    std_dev = np/2;

    % default for std_dev = 128
    a = -log(13.999/1023)/std_dev^2; % to get the exact width of the gaussian in Variam archives
    for j = 1:np
        amplitude = 1023*exp(-a*x^2);
        shape_gauss(j,1) = phase;
        shape_gauss(j,2) = amplitude;
        shape_gauss(j,3) = time_part;
        x = x + dx;
    end
    
    pulse_shape = shape_gauss;

end

function pulse_shape = make_FourierProfile(A_Coef, P_Coef, np)
    % give the pulse profile as a matrix [amplitude phase time_partition] <<<<<
    % A/P_Coef is the Nx3 coefficient such that C(k) = (a_k,b_k,c_k), where 
    %       Omega(t) = \sum(a_k*sin(b_k*t + c_k)) 
    % A_Coef and P_Coef have size s_a and s_p respectively 
    % np is the number of points calculated, i.e, len(time_partition)
    % Modulation inspired by [Peterson et al. Phys. Rev. Applied 13, 054060 (2020)]
   
    t_0 = 0;
    
    Amp = zeros(np,1); Phase = zeros(np,1); time_part = ones(np,1); 

    t = t_0; dt=4.0*pi/np; dt = 1.0;
    for j = 1:np
        for k = 1:size(A_Coef,1)
            Amp(j) = Amp(j) + A_Coef(k,1)*sin(A_Coef(k,2)*t + A_Coef(k,3));
        end
        for k = 1:size(P_Coef,1)
            Phase(j) = Phase(j) + P_Coef(k,1)*sin(P_Coef(k,2)*t + P_Coef(k,3));
        end    
        t = t + dt; 
    end
    tau_f =  t - dt;

    
    
    % Amplitude must be bounded by [0,MaxAmp of our machine]
    % for MaxAmp do: if Amp_max/MaxAmp>1 --> Amp = Amp*MaxAmp/Amp_max; 
    % else, Amp=Amp
    Amp_min = min(Amp);
    if Amp_min < 0
        Amp = Amp - Amp_min*ones(size(Amp,1),1);
    end
    MaxAmp = 1023; 
    Amp = Amp*MaxAmp;
    ratio = max(Amp)/MaxAmp;
    if ratio > 1
        Amp = Amp*MaxAmp/max(Amp);
    end

    % Also Phase must be bounded by [0,2*Pi)
    if Phase > 360
        Phase = Phase/max(Phase)*360;
    end
    Phase = mod(Phase, 360);
    
    % Calculating the pulse envelope:
    zeta_1 =2; zeta_2 = 2; % < ---- these are the Laflamme paper values, others could work also!
    t = t_0; 

    for j = 1:np
        Delta(j) = -tanh(zeta_1*t/tau_f)*tanh(zeta_2*(t - tau_f)/tau_f);
        Amp(j) = Amp(j)*Delta(j);
        t = t + dt;
    end
    Amp = Amp/max(Delta);
    
    pulse_shape = [Phase Amp time_part];

end 

function result = find_FSModulation_loss1q_Hom(U_target,Coeff,s_a,s_p,np, power, pulse_width, rho_ini, L_Spec, L_Mol, isPar)  
    % Given a pulse modulation, this function finds the best Coeff of a Fourier Series so it equals a U_target.
    % --- this is a function for homonuclear systems ---- %
    % L_Spec, L_Mol are matrix versions of global variables (structures) mol and spectro, used only when parallelization is on
    arguments
        U_target (:,:) % unitary evolution target
        Coeff (:,:) % cat(1, Amplitude_Coefficients, Phase_Coefficients) 
        % is matrix Amplitude and Phase congatenated , with size ((s_a+s_p)*Nq,3)
        s_a double = 1 % number of coeff to sum in amplitude
        s_p double = 1 % number of coeff to sum in phase
        np double =  100 % number of points
        power double = 20*1e3 % max_power on channel
        pulse_width = double % time duration of pulse, in (us) 
        rho_ini (:,:) = [1 0; 0 0] % inicial state
        L_Spec (1,:) =  [10.0000   24.8750] 
        L_Mol (:,:,:) = zeros(8,8,6)
        isPar double = 0 %set to 0 when not parallelizing, and 1 otherwise
    end

    % Amp and Phase Coefficients spliting 
    A_Coef = Coeff(1:s_a, :);
    P_Coef = Coeff(s_a+1:s_a + s_p, :);

    shape_mod = make_FourierProfile(A_Coef, P_Coef, np); % generate shape by fourier series
    
    % apply generated shape giving the corresponding unitary and final state
    %[U rohfin] = shapedpulse2(shape,pw,phase,ramp,power,roha);
    if isPar== 0
        [U, rhofin] = shapedpulse_Hom_Par(shape_mod,pulse_width,shape_mod(:,1),0,power,rho_ini, L_Spec, L_Mol); 
    elseif isPar== 1
        [U, rhofin] = shapedpulse_Hom_Par(shape_mod,pulse_width,shape_mod(:,1),0,power,rho_ini, L_Spec, L_Mol); 
    end
    % Function to be minimized chi^2 
    %chi_2 = sum(sum(abs(abs(U) - abs(U_target))));
    chi_2 = sum(sum( real(U - U_target).^2 + imag(U - U_target).^2 ));
    
    result = chi_2;
end

function result = find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pulse_width, ...
                                                    rho_ini, L_Spec, L_Mol, isPar)  
    % Given a Amp pulse modulation, this function finds the best Coeff of a Fourier Series so it equals a U_target.
    % --- this is a function for homonuclear systems ---- %
    % L_Spec, L_Mol are matrix versions of global variables (structures) mol and spectro, used only when parallelization is on
    arguments
        U_target (:,:) % unitary evolution target
        Coeff (:,:) % cat(1, Amplitude_Coefficients, Phase_Coefficients) 
        % is matrix Amplitude and Phase congatenated , with size ((s_a+s_p)*Nq,3)
        AmpShape (:,3)  % amplitude as a matrix, usually taken from an .RF file
        s_p double = 1 % number of coeff to sum in phase
        np double =  100 % number of points
        power double = 20*1e3 % max_power on channel
        pulse_width = double % time duration of pulse, in (us) 
        rho_ini (:,:) = [1 0; 0 0] % inicial state
        L_Spec (1,:) =  [10.0000   24.8750] 
        L_Mol (:,:,:) = zeros(8,8,6)
        isPar double = 0
    end
    
    % Phase Coefficients 
    P_Coef = Coeff(1:s_p, :);
    A_Coef = zeros(1,3);

    shape_mod = make_FourierProfile(A_Coef, P_Coef, np); % generate shape by fourier series
    
    % apply generated shape giving the corresponding unitary and final state
        %[U, rhofin] = shapedpulse_Hom(AmpShape,pulse_width,shape_mod(:,1),0,power,rho_ini); %, L_Spec, L_Mol); 
    if isPar == 1
        [U, rhofin] = shapedpulse_Hom_Par(AmpShape,pulse_width,shape_mod(:,1),0,power,rho_ini, L_Spec, L_Mol); 
    else
        [U, rhofin] = shapedpulse_Hom(AmpShape,pulse_width,shape_mod(:,1),0,power,rho_ini); %, L_Spec, L_Mol); 
    end

    % Function to be minimized chi^2 
    %chi_2 = sum(sum(abs(abs(U) - abs(U_target))));
    chi_2 = sum(sum( real(U - U_target).^2 + imag(U - U_target).^2 ));
    
    result = chi_2;
end

function result = find_phase_FSModulation_loss1q_Het(chn, U_target,Coeff,AmpShape_obs, AmpShape_dec,s_p,np, ...
                                power_obs, power_dec, pw_obs, pw_dec, rho_ini, L_Spec, L_Mol, isPar)  
    % Given a Amp pulse modulation, this function finds the best Coeff of a Fourier Series so it equals a U_target.
    % --- this is a function for heteronuclear systems ---- %
    % L_Spec, L_Mol are matrix versions of global variables (structures) mol and spectro, used only when parallelization is on
    arguments
        chn double % either 0 or 1 for obs/dec
        U_target (:,:)  % unitary evolution target
        Coeff (:,:) % , Phase_Coefficients) 
        % is matrix Amplitude and Phase congatenated , with size ((s_a+s_p)*Nq,3)
        AmpShape_obs (:,3)  % amplitude as a matrix, usually taken from an .RF file
        AmpShape_dec (:,3)  % amplitude as a matrix, usually taken from an .RF file
        s_p double = 1 % number of coeff to sum in phase
        np double =  100 % number of points
        power_obs double = 20*1e3 % max_power on channel obs
        power_dec double = 20*1e3 % max_power on channel dec
        pw_obs double = 200 % time duration of pulse, in (us) 
        pw_dec double = 200
        rho_ini (:,:) = [1 0; 0 0] % inicial state
        L_Spec (1,:) =  [10.0000   24.8750] 
        L_Mol (:,:,:) = zeros(8,8,6)
        isPar double = 0
    end
    
    % Phase Coefficients 
    P_Coef = Coeff(1:s_p, :);
    A_Coef =[];

    shape_mod = make_FourierProfile(A_Coef, P_Coef, np); % generate shape by fourier series
    AmpShape_obs(:,1) = shape_mod(:,1);
    AmpShape_dec(:,1) = shape_mod(:,1);
    

    % apply generated shape giving the corresponding unitary and final state
        %[U, rhofin] = shapedpulse_Hom(AmpShape,pulse_width,shape_mod(:,1),0,power,rho_ini); %, L_Spec, L_Mol); 
    if isPar == 1
        %[U, roh] = shapedpulse_Par(shape_obs,pw_obs,phase_obs,ramp_obs,power_obs,roha, L_Spec, L_Mol,channel)
        if chn == 0
            [U, rhofin] = shapedpulse_Par(AmpShape_obs, pw_obs, [], 0, power_obs,rho_ini, L_Spec, L_Mol, chn); 
        elseif chn ==1
            [U, rhofin] = shapedpulse_Par(AmpShape_dec, pw_dec, [], 0, power_dec,rho_ini, L_Spec, L_Mol, chn); 
        end
    else
        if chn == 0
            [U, rhofin] = shapedpulse_Par(AmpShape_obs, pw_obs, [], 0, power_obs,rho_ini, L_Spec, L_Mol, chn); 
        elseif chn ==1
            [U, rhofin] = shapedpulse_Par(AmpShape_dec, pw_dec, [], 0, power_dec,rho_ini, L_Spec, L_Mol, chn); 
        end
    end

    % Function to be minimized chi^2 
    chi_2 = sum(sum( real(U - U_target).^2 + imag(U - U_target).^2 ));
    
    result = chi_2;
end

function result = find_phase_FSModulation_loss2q_Het(U_target,Coeff,AmpShape_obs, AmpShape_dec,s_p_obs, s_p_dec,np, ...
                                power_obs, power_dec, pw_obs, pw_dec, rho_ini, L_Spec, L_Mol, isPar)  
    % Given a Amp pulse modulation, this function finds the best Coeff of a Fourier Series so it equals a U_target.
    % --- this is a function for heteronuclear systems ---- %
    % L_Spec, L_Mol are matrix versions of global variables (structures) mol and spectro, used only when parallelization is on
    arguments
        U_target (:,:) % unitary evolution target
        Coeff (:,:) % cat(1, Amplitude_Coefficients, Phase_Coefficients) 
        % is matrix Amplitude and Phase congatenated , with size ((s_a+s_p)*Nq,3)
        AmpShape_obs (:,3)  % amplitude as a matrix, usually taken from an .RF file
        AmpShape_dec (:,3)  % amplitude as a matrix, usually taken from an .RF file
        s_p_obs double = 1 % number of coeff to sum in phase
        s_p_dec double = 1 % number of coeff to sum in phase
        np double =  100 % number of points
        power_obs double = 20*1e3 % max_power on channel obs
        power_dec double = 20*1e3 % max_power on channel dec
        pw_obs double = 200% time duration of pulse, in (us) 
        pw_dec double = 200
        rho_ini (:,:) = [1 0; 0 0] % inicial state
        L_Spec (1,:) =  [10.0000   24.8750] 
        L_Mol (:,:,:) = zeros(8,8,6)
        isPar double = 0
    end
    
    % Phase Coefficients 
    P_Coef_obs = Coeff(1:s_p_obs, :);
    P_Coef_dec = Coeff(s_p_obs+ 1: s_p_obs + s_p_dec, :);
    A_Coef =[];

    shape_mod_obs = make_FourierProfile(A_Coef, P_Coef_obs, np); % generate shape by fourier series
    shape_mod_dec = make_FourierProfile(A_Coef, P_Coef_dec, np); % generate shape by fourier series
    AmpShape_obs(:,1) = shape_mod_obs(:,1);
    AmpShape_dec(:,1) = shape_mod_dec(:,1);

    % apply generated shape giving the corresponding unitary and final state
        %[U, rhofin] = shapedpulse_Hom(AmpShape,pulse_width,shape_mod(:,1),0,power,rho_ini); %, L_Spec, L_Mol); 
    if isPar == 1
        %[U, roh] = simshapedpulse2_Par(shape_obs,shape_dec,pw_obs,pw_dec,... 
        %            ramp_obs,ramp_dec,power_obs,power_dec,roha, L_Spec, L_Mol)
        [U, rhofin] = simshapedpulse2_Par(AmpShape_obs, AmpShape_dec,pw_obs, pw_dec, 0, 0, ...
                            power_obs,power_dec,rho_ini, L_Spec, L_Mol); 
    else
        [U, rhofin] = simshapedpulse2_noPar(AmpShape_obs, AmpShape_dec,pw_obs, pw_dec, 0, 0, ...
                            power_obs,power_dec,rho_ini);
    end

    % Function to be minimized chi^2 
    chi_2 = sum(sum( real(U - U_target).^2 + imag(U - U_target).^2 ));
    
    result = chi_2;
end

function result = find_FSModulation_loss1q_Het(...
        U_target,Coeff,s_a,s_p,np, power, pulse_width, rho_ini ,L_Spec, L_Mol, isPar,channel) 
    % Given a pulse modulation, this function finds the best Coeff of a Fourier Series so it equals a U_target.
    % --- this is a function for heteronuclear systems target a frequency-lone qbit ---- %
    %{
    arguments
        U_target (:,:) % unitary evolution target
        Coeff (:,:) % cat(1, Amplitude_Coefficients, Phase_Coefficients) 
        % is matrix Amplitude and Phase congatenated , with size ((s_a+s_p)*Nq,3)
        s_a double = 1 % number of coeff to sum in amplitude
        s_p double = 1 % number of coeff to sum in phase
        np double =  100 % number of points
        power double = 20*1e3 % max_power on channel
        pulse_width = double % time duration of pulse, in (us) 
        rho_ini (:,:) = [1 0; 0 0] % inicial state
    end
    %}
    % Amp and Phase Coefficients spliting 
    A_Coef = Coeff(1:s_a, :);
    P_Coef = Coeff(s_a+1:s_a + s_p, :);
    
    % generate modulated pulse by fourier series
    shape_mod = make_FourierProfile(A_Coef, P_Coef, np);
    
    % generate the corresponding unitary and final state
    %[U roh] = shapedpulse2(shape,pw,phase,ramp,power,roha);
    [U, rhofin] = shapedpulse_Par(shape_mod,pulse_width,[],0,power,rho_ini, L_Spec,L_Mol,channel);
    
    % Function to be minimized chi^2 
    chi_2 = sum(sum( real(U - U_target).^2 + imag(U - U_target).^2 ));
    
    result = chi_2;

end

function result = find_FSModulation_loss2q(U_target,Coeff_ini2q,s_aVec,s_pVec,np, power_obs,power_dec, pulse_width, rho_ini, L_Spec, L_Mol, isPar) 
    % Given a pulse modulation, this function finds the best Coeff of a Fourier Series so it equals a U_target.
    % --- this is a function for heteronuclear systems target two separeted frequency qbits ---- %
    arguments
        U_target (:,:) % unitary evolution target
        Coeff_ini2q (:,:) % cat(1, Amplitude_Coefficients, Phase_Coefficients) 
        % is matrix Amplitude and Phase congatenated , with size ((s_a+s_p)*Nq,3), Nq= number of qbits
        s_aVec (1,:) % number of coeff to sum in amplitude
        s_pVec (1,:) % number of coeff to sum in phase
        np double =  100 % number of points
        power_obs double = 20*1e3 % max_power on observer channel 
        power_dec double = 20*1e3 % max_power on observer channel 
        pulse_width = double % time duration of pulse, in (us) 
        rho_ini (:,:) = [1 0; 0 0] % inicial state
        L_Spec (1,:) =  [10.0000   24.8750] 
        L_Mol (:,:,:) = zeros(8,8,6)
        isPar double = 0
    end
    s_p_obs = s_pVec(1);
    s_p_dec = s_pVec(2);
    s_a_obs = s_aVec(1);
    s_a_dec = s_aVec(2);
    sum_sa = s_a_obs +s_a_dec;
    sum_sp = s_p_obs+ s_p_dec; 

    % Amp and Phase Coefficients spliting 
    A_Coef_obs = Coeff_ini2q(1:s_a_obs, :);
    P_Coef_obs = Coeff_ini2q(s_a_obs+1:s_a_obs + s_p_obs, :);
    A_Coef_dec = Coeff_ini2q(s_a_obs + s_p_obs + 1 : sum_sa + s_p_obs, :);
    P_Coef_dec = Coeff_ini2q(sum_sa + s_p_obs + 1 : sum_sa +sum_sp, :);

    % generate modulated pulse by fourier series
    shape_mod_obs = make_FourierProfile(A_Coef_obs, P_Coef_obs, np);
    shape_mod_dec = make_FourierProfile(A_Coef_dec, P_Coef_dec, np);
    

    % generate the corresponding unitary and final state
    %[U, roh] = simshapedpulse2_Par(shape_obs,shape_dec,pw_obs,pw_dec,... 
    %ramp_obs,ramp_dec,power_obs,power_dec,roha, L_Spec, L_Mol)
    if isPar == 1
        [U, roh] = simshapedpulse2_Par(shape_mod_obs,shape_mod_dec,pulse_width,pulse_width,...
            0,0,power_obs,power_dec,rho_ini, L_Spec, L_Mol);
    elseif isPar ==0
        [U, roh] = simshapedpulse2_noPar(shape_mod_obs,shape_mod_dec,pulse_width,pulse_width,...
            0,0,power_obs,power_dec,rho_ini); %, L_Spec, L_Mol);
    end

    % Function to be minimized chi^2 
    chi_2 = chi2_loss(U,U_target);
    %noise = [0.05 0.05 0.05 0.05 0.05]; % hardcoded noise
    %result = noisy_loss(U,U_target,rho_ini,noise);
    %theta = angle(U_target(1,1)) - angle(U(1,1));
    %chi_2 = sum(sum( real(U - U_target).^2 + imag(U - U_target).^2 )) + abs(theta);
    result = chi_2;
    
end


function chi_2 = chi2_loss(U_gen, U_target)
    % loss function for noiseless qubits

    theta = angle(U_target(1,1)) - angle(U_gen(1,1));
    % Function to be minimized chi^2 
    chi_2 = sum(sum( real(U_gen - U_target).^2 + imag(U_gen - U_target).^2 )) + abs(theta);
end

function loss = noisy_loss(U_gen,U_target, rhoini, noise)
    %calculates the loss given a U_target and expected noise
    % multiplies the states as:
    % finds number of qbits
    [a1,~] = size(rhoini);
    a11 = log2(a1);
    input_states = multiplicate_states(rhoini,a11);
    p = noise(1:4);
    q = noise(5);
    chi_2 = 0;
    for j=1:a11
        for i=1:length(input_states)
        output_states_target{j}{i} = PartialTrace(ApplyOperator(input_states{i}, U_target),j);
        GenOutput_states_target{j}{i} = PartialTrace(ApplyOperator(input_states{i}, U_gen),j);
        GenOutput_states_target{j}{i} = Apply_composite_map(GenOutput_states_target{j}{i}, p, q);
        Delta = output_states_target{j}{i} - GenOutput_states_target{j}{i};
        loss = sum(sum( real(Delta).^2 + imag(Delta).^2 ));
        end
    end
end

function products = multiplicate_states(rho,numSubsystems)
    % Define Pauli matrices
    Pauli1 = [0 1; 1 0];
    Pauli2 = [0 -1i; 1i 0];
    Pauli3 = [1 0; 0 -1];
    Identity = eye(2);
    Paulis = {Identity, Pauli1, Pauli2, Pauli3};

    % Generate all possible combinations of Pauli matrices for the given number of subsystems
    [indices{1:numSubsystems}] = ndgrid(1:4);
    indices = cellfun(@(x) x(:), indices, 'UniformOutput', false);
    allCombinations = [indices{:}];

    % Compute product for each combination and store in products cell array
    products = cell(size(allCombinations, 1), 1);
    for i = 1:size(allCombinations, 1)
        combination = allCombinations(i, :);
        operator = Paulis{combination(1)};
        for j = 2:numSubsystems
            operator = kron(operator, Paulis{combination(j)});
        end
        products{i} = operator * rho * operator';
    end
end

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

function loss = offset_lossHom(PulseMod, gate_num, offset, rho_ini, L_Spec) 
        currentQbit = PulseMod.QbitOrder{gate_num};
        U_target = PulseMod.U_target{gate_num};
        pw_obs = PulseMod.width_obs{gate_num};
        power_bound_obs = PulseMod.power_bound_obs{gate_num};
        power = power_bound_obs(2);
        s_a = PulseMod.s_a_obs{gate_num};
        s_p = PulseMod.s_p_obs{gate_num};
        np = PulseMod.np_obs{gate_num};
        Par = PulseMod.Par;
        isPar = Par(1);
        AmpShape = PulseMod.shape_obs{gate_num};

        L_Mol = focusfrequency(currentQbit, offset);
        loss = find_phase_FSModulation_loss1q_Hom(U_target,Coeff,AmpShape,s_p,np, power, pw_obs, ...
                                                            rho_ini, L_Spec, L_Mol, isPar)  ;
        L_Mol = focusfrequency(target_qbit, offset=0);

end

