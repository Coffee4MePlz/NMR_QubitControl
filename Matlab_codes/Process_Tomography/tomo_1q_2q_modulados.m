% -------------------------------------------------------------------------
% Tomografia 2  qubits 2 qubit - JOB:  Pulsos Modulados cafe
% May 26, 2023
% -------------------------------------------------------------------------
%
clc
clear all
close all

% Start - adiciona os paths do scripts importantes
%{
addpath(['/home/carlos/MEGAsync/NMR_NEW/Matlab_codes/nmr_lib/'])    % scripts para simulacao da evolucao e pulsos NMR
addpath(['/home/carlos/MEGAsync/NMR_NEW/Matlab_codes/gates/'])      % gates
addpath(['/home/carlos/MEGAsync/NMR_NEW/Matlab_codes/Math_util/'])  % math utils 
addpath(['/home/carlos/MEGAsync/NMR_NEW/Matlab_codes/molecules/'])  % dados das moleculas
addpath(['/home/carlos/MEGAsync/NMR_NEW/Matlab_codes/control/'])
addpath(['/home/carlos/MEGAsync/NMR_NEW/Matlab_codes/Carlos/'])
actdir =  pwd;
%}
%%

% Definicoes uteis
id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
ZZ = kron(sz,sz);

% funções de rotação angulo thetai
Rotx = @(thet1) [cos(thet1/2.0) -1i*sin(thet1/2.0);  -1i*sin(thet1/2.0) cos(thet1/2.0)];
Roty = @(thet2) [cos(thet2/2.0) -sin(thet2/2.0);  sin(thet2/2.0) cos(thet2/2.0)];
Rotz = @(thet3) [exp(-1i*thet3/2.0) 0; 0 exp(1i*thet3/2.0)];
% projetores z
prj{1} = [1.0 0.0; 0.0 0.0];
prj{2} = [0.0 0.0; 0.0 1.0];
% projetores x
prjx{1} = 0.5*[1.0 1.0; 1.0 1.0];      %|+><+|
prjx{2} = 0.5*[1.0 -1.0; -1.0 1.0];    %|-><-|
% projetores y
prjy{2} = 0.5*[1.0 -1i; 1i 1.0];      %|+><+|
prjy{1} = 0.5*[1.0 1i; -1i 1.0];    %|-><-|

% -------
spectropar; % Parametros do espectometro
cloroformio; % Parametrso da molecula - formato de sódio - H C O_2 Na

% ---------------------------------------------------------------------------
% >>>>>>>>>>>>    Primeiro passo - Espectro do equilibrio nr = 0
% Ajustar os parametros do fit
% Encontra o eixo da frequencia (ff1 - ressonancia)
[tran1 tran1 ff1 , opr] = transitions(mol,1);
[tran2 tran2 ff2 opr] = transitions(mol,2);

% Medindo H
ff1 = ff1-0.054; %para testar se o primeiro qubit esta fora da ressonancia resolucao espectometro +- 0.10 Hz
t2 = 0.617;     %t2 para tetar o fit 

% Atualizar o nome do arquivo fid sem eqtest.fido ponto
% Carbono tipicamente 64*1024 hidrogenio 32*1024
spechr= varian2mat('exp400',32*1024); % feito no dia 27/05-23

% Fita e plota o espectro - mexer no t2 para ajustar o fit
% fitespec(spechr{1},ff1,t2,[intervalo para freq. ressonancia, intervalo para t2],[comprimento do eixo x - freuqencia relativa]);
spechr{1} = fitespec(spechr{1},ff1,t2,[0.005 0.005],[-200 200]);
plotfit(1,spechr{1},[-200 200]); 

% Calc. acoplamento J
Jint_from_spec = spechr{1,1}.fr(1)-spechr{1,1}.fr(2);

%% Para plotar o fid
%t2_from_fid = plota_fid_T2(spechr,2)

%% ---------------------------------------------------------------------------
% >>>>>>>>>>>> Segundo passo - Inicio das tomografias de 2 qubits
%% ---------------------------------------------------------------------------
% Tomo do |00> - calculo da normalizaçào padrão para o PPS

aa=varian2mat('exp405',64*1024); % Parte 1 da tomo
bb=varian2mat('exp406',64*1024); % Parte 2 da tomo

[rho00, desv_00, norma2q] = tomografia_00(aa,bb,ff1,t2,0);

target = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

PrettyPlotMatrix3(rho00{1},30,1)

% Autovalores de roh
autoval00 = eig(rho00{1});
% Matriz Corrigida Max. Likelly.
rho_ml_00 = find_rho_2q(rho00{1},0);
% Calcula a fidelidade
fidd_00 =  fidelmat(rho00{1},target);
% Distancia de traco
tracedist_00 = trdis(rho00{1},target);
% Calcula a fidelidade rohml
fiddml_00 =  fidelmat(rho_ml_00,target);
% Distancia de traco rohml
tracedistml_00 = trdis(rho_ml_00,target);
% Uhlmann fidelity of two density matrices
Uhlmann_fidelity_00 = Fidelity(target,rho_ml_00);

%% ---------------------------------------------------------------------------
% Tomo de 1 qubit do |00> - calculo da normalizaçào padrão para o PPS

aa_aux=varian2mat('exp405',64*1024); % Parte 1 da tomo
aa_1q{1} = aa_aux{1};
aa_1q{2} = aa_aux{2};

[rho_1q_00, desv_1q_00, norma1q] = tomografia_1q_00(aa_1q,ff1,t2,0);

target = [1 0; 0 0];

showmat(33,rho_1q_00{1});

%PrettyPlotMatrix3(rho00{1},30,1)

% Autovalores de roh
autoval_1q_00 = eig(rho_1q_00{1});
% Calcula a fidelidade
fidd_1q_00 =  fidelmat(rho_1q_00{1},target);


%% ---------------------------------------------------------------------------
% Tomo de 1 qubit da superposicao  ?

aa_1q=varian2mat('exp407',32*1024); % Parte 1 da tomo Rota 1
[rho_1q_10, desv_1q_10] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);
% Para normalizar o |1> 
norma1q_10 = (1/2)/desv_1q_10{1}(2,2);
% Estado normalizado
rho_1q_10{1} = (norma1q_10).*desv_1q_10{1} + eye(2)/2;

%%
aa_1q=varian2mat('exp408',32*1024); % Parte 1 da tomo Rota 1
[rho_1q_xplus, desv_1q_xplus] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);
rho_1q_xplus{1}
auxx = Roty(pi/2)*[1,0;0,0]*Roty(-pi/2);
fidelmat(rho_1q_xplus{1},auxx)
%trace distance
trdis(rho_1q_xplus{1},auxx)
% Uhlmann fidelity of two density matrices
Fidelity(rho_1q_xplus{1},auxx)
%%
aa_1q=varian2mat('exp409',64*1024); % Parte 1 da tomo Rota 1
[rho_1q_yminus, desv_1q_yminus] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);
rho_1q_yminus{1}
%% Estados apos o pulso pi/2
aa_1q=varian2mat('exp410',64*1024); % Parte 1 da tomo Rota 1
[rho_1q_Rot00, desv_1q_Rot00] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);
rho_1q_Rot00{1}
%%
aa_1q=varian2mat('exp411',64*1024); % Parte 1 da tomo Rota 1
[rho_1q_Rot10, desv_1q_Rot10] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q_10,0);
rho_1q_Rot10{1}
%% 
aa_1q=varian2mat('exp412',64*1024); % Parte 1 da tomo Rota 1
[rho_1q_Rotxplus, desv_1q_Rotxplus] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);

rho_1q_Rotxplus{1}

%% 
aa_1q=varian2mat('exp413',64*1024); % Parte 1 da tomo Rota 1
[rho_1q_Rotyminus, desv_1q_Rotyminus] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);

rho_1q_Rotyminus{1}


%%
save(['/home/carlos/MEGAsync/NMR_NEW/Matlab_codes/Carlos/rho_tomo_proc'],'rho_1q_00','rho_1q_10','rho_1q_xplus',...
    'rho_1q_yminus','rho_1q_Rot00','rho_1q_Rot10','rho_1q_Rotxplus','rho_1q_Rotyminus');


%% 
aa_1q=varian2mat('exp420',32*1024); % Parte 1 da tomo Rota 1
[rho_1q_Rotxmod, desv_1q_Rotxmod] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);

rho_1q_Rotxmod{1}

%%
aa_1q=varian2mat('exp432',32*1024); % Parte 1 da tomo Rota 1
[rho_1q_Rotymod, desv_1q_Rotymod] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);

rho_1q_Rotymod{1}
fidelmat(rho_1q_Rotymod{1}, prjy{2})

%% 
% Fullmodulados -- simultâneo Roty(90)
TomoStates{1} = prjx{1};
TomoStates{2} = prjx{2};
TomoStates{3} = prj{2};
TomoStates{4} = prjy{2};
F_mod2 = 1;
for j=1:4
aa_1q=varian2mat(append('exp41', num2str(j-1)),32*1024); % Parte 1 da tomo Rota 1
[rho_1q_Rotymod2{j}, desv_1q_Rotymod2] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);
    F_mod2(j) = fidelmat(rho_1q_Rotymod2{j}{1}, TomoStates{j}) ;
end

F_mod2
rho_1q_Rotymod2{1}{1}

%auxx = Roty(pi/2)*[1,0;0,0]*Roty(-pi/2);
%fid_full = fidelmat(rho_1q_Rotymod2{1},auxx)
%trace distance
%trdis(rho_1q_Rotymod2{1},auxx)
% Uhlmann fidelity of two density matrices
%Fidelity(rho_1q_Rotymod2{1},auxx)

%%
% ----- Tomo do CNOT|00> - pulso hard ----- %%%

aa=varian2mat('exp476',64*1024); % Parte 1 da tomo
bb=varian2mat('exp477',64*1024); % Parte 2 da tomo

[rhocnot, desv_00, norma2q] = tomografia_00(aa,bb,ff1,t2,0);
%[rhocnot, desv_00] = tomografia_1q_mult(aa,ff1,t2,norma1q,0);

%target2 = [1 0 ; 0 0];
target2 = Cnot1([1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]);

PrettyPlotMatrix3(rhocnot{1},50,1)

fidelmat(rhocnot{1},target2)

%%
%%
% ----- Tomo do CNOT|00> - modulado ----- %%%


aa=varian2mat('exp480',32*1024); % Parte 1 da tomo
bb=varian2mat('exp481',32*1024); % Parte 2 da tomo

[rhocnot, desv_00, norma2q] = tomografia_00(aa,bb,ff1,t2,0);

[rhocnot_1q, desv_00] = tomografia_1q_mult(aa,ff1,t2,norma1q,0);

target2 = Cnot1([1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]);
target1 = [1 0 ; 0 0];


PrettyPlotMatrix3(rhocnot{1},50,1)
PrettyPlotMatrix3(rhocnot_1q{1},51,1)

fidelmat(rhocnot{1},target2)
fidelmat(rhocnot_1q{1},target1)

%%
% ----- Tomo do CNOT|00> - modulado ----- %%%

%TomoStates2q{1} = kron(prj{1},prj{1});
%TomoStates2q{2} = kron(prj{2},prj{2});
TomoStates2q{1} = CNOT*kron(prjx{1},prj{1})*CNOT';
TomoStates2q{2} = CNOT*kron(prjy{2},prj{1})*CNOT';
%TomoStates2q{5} = kron(prj{2},prj{2});

%TomoStates1q{1} = prj{1};
%TomoStates1q{2} = prj{2};
TomoStates1q{1} = prjx{1};
TomoStates1q{2} = prjy{2};
%TomoStates1q{5} = prj{2};

F_mod_1q =1;
F_mod =1;
for j=1:2
aa=varian2mat(append('exp49', num2str(2*(j-1))),32*1024); % Parte 1 da tomo
bb=varian2mat(append('exp49', num2str(2*j -1)),32*1024); % Parte 2 da tomo

[rhocnot{j}, desv_00, norma2q] = tomografia_00(aa,bb,ff1,t2,0);
[rhocnot_1q{j}, desv_00] = tomografia_1q_mult(aa,ff1,t2,norma1q,0);
F_mod(j) = fidelmat(rhocnot{j}{1},TomoStates2q{j});
F_mod_1q(j) = fidelmat(rhocnot_1q{j}{1},TomoStates1q{j});
end
%target2 = [1 0 ; 0 0];

target2 = Cnot1([1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]);
target1 = [1 0 ; 0 0];


PrettyPlotMatrix3(rhocnot{1}{1},50,1)
PrettyPlotMatrix3(rhocnot_1q{1}{1},51,1)

fidelmat(rhocnot{1}{1},target2)
fidelmat(rhocnot_1q{1}{1},target1)
%%

%%
% Fullmodulados X pulso hard // teste para Offset
aa_1q=varian2mat('exp500',32*1024); 
[rho_1q_Rotyhard2, desv_1q_Rotyhard2] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);

aa_1q=varian2mat('exp501',32*1024); 
[rho_1q_Rotymod2, desv_1q_Rotymod2] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);

rho_1q_Rotyhard2{1};

auxx = Roty(pi/2)*[1,0;0,0]*Roty(-pi/2);
for j=1:21
    disp('hard pulse w/ offset:')
    F_hard(j) = fidelmat(rho_1q_Rotyhard2{j},auxx);
    disp('mod pulse w/ offset')
    F_mod(j) = fidelmat(rho_1q_Rotymod2{j}, auxx);
end

    disp('hard pulse w/ offset:')
    F_hard
    disp('mod pulse w/ offset')
    F_mod
%% testes de pulsos a tempos maiores:
% exp502-7, onde impares são Roty(90) e pares são Rotx(90)

for j=1:6
    str = append('exp50',num2str(j+1));
aa_1q=varian2mat(str,32*1024); 
[rho_1q_Rotxymod{j}, desv_1q_Rotxymod{j}] = tomografia_1q_mult(aa_1q,ff1,t2,norma1q,0);
end
%%
targety = expm(-1i*pi/4*sx)*prj{1}*expm(1i*pi/4*sx);
targetx = expm(-1i*pi/4*sy)*prj{1}*expm(1i*pi/4*sy);
Fid_200 = [fidelmat(rho_1q_Rotxymod{1}{1}, targety), fidelmat(rho_1q_Rotxymod{2}{1}, targety), fidelmat(rho_1q_Rotxymod{3}{1}, targety)]
Fid_400 = [fidelmat(rho_1q_Rotxymod{4}{1}, targety), fidelmat(rho_1q_Rotxymod{5}{1}, targety), fidelmat(rho_1q_Rotxymod{6}{1}, targety)]

%%

% ------------------------------------------------------------------------
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ------------------------------------------------------------------------
% CNOT - controle qubit 1
function result = Cnot1(rho)
id = [1.0 0.0; 0.0 1.0];
sx = [0.0 1.0; 1.0 0.0];
% Cnot = |0><0|x I + |1><1|x sigma_x
aux = kron([1.0 0.0; 0.0 0.0],id) + kron([0.0 0.0; 0.0 1.0],sx);
result = aux*rho*ctranspose(aux);
end
% ---------------------------------------------------------------------------
function rho_out = gradiente(rho)
% Simula o gradiante de forma simplificada
% ---------------------------------------------------------------------------

% Gradiente 
for i = 1:4
  for j = 1:4
    if i ~= j
      rho(i,j) = rho(i,j)/20;
      %rho(i,j) = 0.0;
    end
  end
end

rho_out = rho;

end

% ---------------------------------------------------------------------------
function t2_from_fid = plota_fid_T2(spechr,nfig)
% Plota o fid do espectro de equilibrio e determinando t2 a partir do fid
% ---------------------------------------------------------------------------

signal_fid = real(spechr{1,1}.fid); % Parte real do Fid
atemp = linspace(0,8,length(signal_fid)); %  preparando o eixo do tempo 0 a 4s ou 0 a 8s

figure(nfig)
plot(atemp,signal_fid) % plota o fid em unidades arbitrarias
figpar = gca(figure(21));
figpar.FontSize = nfig;
ylabel('FID (a.u.)','fontsize',18)
xlabel('t ($\mu$s)',...
    'Interpreter','latex','FontSize',18)
xlim([0 4])

% Pegando os picos do fid
[peakValues, peakIndexes] = findpeaks(signal_fid);

%Limpando eventuais picos negativos para plotar o decaimento do fid
jaux=1;
for j = 1:length(peakIndexes)
    if peakValues(j) >= 0
        ay(jaux) = peakValues(j);
        ax(jaux) = atemp(peakIndexes(j));
        jaux=jaux+1;
    end
end

nay = ay/max(ay); % decaimento do fid normalizado

fit_fid_decay=fit(ax',nay','exp2'); % decaimento do fid fitado com uma dupla exponecial

figure(nfig+1)
plot(fit_fid_decay,ax,nay);
figpar = gca(figure(nfig+1));
figpar.FontSize = 14;
ylabel('FID decay','fontsize',18)
xlabel('t (s)',...
    'Interpreter','latex','FontSize',18)
xlim([0 8])
ylim([-0.01 1.01])

figure(nfig+2)
nsignal_fid = signal_fid/max(signal_fid);  % normaliza o fid
% Para ver o ajuste entre o fit do decaimento e o fid
plot(fit_fid_decay,atemp,nsignal_fid,'-');
figpar = gca(figure(nfig+2));
figpar.FontSize = 14;
ylabel('FID decay','fontsize',18)
xlabel('t (s)',...
    'Interpreter','latex','FontSize',18)
xlim([0 2])
ylim([-1.1 1.1])

% Calculando t2 via decaimento do fid | fit_fid_decay(t) = 1/e = 0.3679 |
t2_from_fid = fzero(@(x) (fit_fid_decay(x) - exp(-1)), [0 4]);

end

% ---------------------------------------------------------------------------
function [rho_out, desv_out, norma2q] = tomografia_00(aa,bb,ff1,t2,flag)
% Razao entre as magnetizacoes
% f =4*nc/nh; % medindo em 2 nucleos corrigir porpocao dos espectros
f = 1; % medino no mesmo nucleo
% fit da parte A
for k=1:length(aa);
    aa{k} =  fitespec(aa{k},ff1,t2,[0.05 0.2],[-200 200]);
    if flag == 1 
        plotfit(k+10,aa{k},[-200 200]); % para ver todos os espectros
    end
end
% fit da parte B
for k=1:length(bb); 	
	 bb{k} =  fitespec(bb{k},ff1,t2,[0.05 0.2],[-200 200]);
     if flag == 1
         plotfit(k+20,bb{k},[-200 200]); % para ver todos os espectros
     end
end

% Inicializa contador para fazer mais uma tomografia
nto = 1; % 

for tt = 1:1:length(aa)     
     k = mod(tt,4);
     if k == 0
         k = 4;
     end
     % atribui a estrtura o valor relacionado a magnetizacao para o H
     struc{1,k} =  f*aa{tt}.ILP;
     % atribui a estrtura o valor relacionado a magnetizacao para o C
     struc{2,k} =  bb{tt}.ILP;
     if k == 4
         % Confiriri os pulsos da sequencia de aquisi??o 
         roh{nto} = tomography({'+II','+YI','+XI','+XX'},struc,2);
         % Para normalizar o |0,0> 
         normalizacao2q = (3.0/4.0)/roh{nto}(1,1);
         % Deviation matrix
         adesv{nto} = roh{nto};
      
         % Estado normalizado
         roh{nto} = (normalizacao2q).*roh{nto} + eye(4)/4;

         nto = nto + 1;
     end
end

rho_out = roh; 
desv_out = adesv; 
norma2q = normalizacao2q;

end

%-------------------- 1 qubit
function [rho_out, desv_out, norma1q] = tomografia_1q_00(aa,ff1,t2,flag)
% fit da parte A
for k=1:length(aa)
    aa{k} =  fitespec(aa{k},ff1,t2,[0.05 0.2],[-200 200]);
    if flag == 1 
        plotfit(k+10,aa{k},[-200 200]); % para ver todos os espectros
    end
end

% Inicializa contador para fazer mais uma tomografia
nto = 1; % 
f = 1;

for tt = 1:1:length(aa)     
     k = mod(tt,2);
     if k == 0
         k = 2;
     end
     % atribui a estrtura o valor relacionado a magnetizacao para o H
     struc{1,k} =  f*aa{tt}.ILP;
%      % atribui a estrtura o valor relacionado a magnetizacao para o C
%      struc{2,k} =  bb{tt}.ILP;
     if k == 2
         % Confiriri os pulsos da sequencia de aquisi??o 
         roh{nto} = tomography({'+II','+YI'},struc,1);
         % Para normalizar o |0> 
         normalizacao1q = (1/2)/roh{nto}(1,1);
         % Deviation matrix
         adesv{nto} = roh{nto};
      
         % Estado normalizado
         roh{nto} = (normalizacao1q).*roh{nto} + eye(2)/2;

         nto = nto + 1;
     end
end

% % Alternativamente
% % Magnetizacao do Hidrogenio 
% mxh = abs(aa{k-1}.ILP(1) + aa{k-1}.ILP(3))/2; % II
% myh = abs(aa{k-1}.ILP(2) + aa{k-1}.ILP(4))/2; % II 
% mzh = abs(aa{k}.ILP(1) + aa{k}.ILP(3))/2;     % YI         
% rohalt{nto} = (mxh*sigmax + myh*sigmay + mzh*sigmaz)

rho_out = roh; 
desv_out = adesv; 
norma1q = normalizacao1q;

end

% ---------------------------------------------------------------------------
function [rho_out, desv_out] = tomografia_multi(aa,bb,ff1,t2,normalizacao2q,flag)
% Razao entre as magnetizacoes
% f =4*nc/nh; % medindo em 2 nucleos corrigir porpocao dos espectros
f = 1; % medino no mesmo nucleo
% fit da parte A
for k=1:length(aa);
    aa{k} =  fitespec(aa{k},ff1,t2,[0.05 0.2],[-200 200]);
    if flag == 1 
        plotfit(k+10,aa{k},[-200 200]); % para ver todos os espectros
    end
end
% fit da parte B
for k=1:length(bb); 	
	 bb{k} =  fitespec(bb{k},ff1,t2,[0.05 0.2],[-200 200]);
     if flag == 1
         plotfit(k+20,bb{k},[-200 200]); % para ver todos os espectros
     end
end

% Inicializa contador para fazer mais uma tomografia
nto = 1; % 

for tt = 1:1:length(aa)     
     k = mod(tt,4);
     if k == 0
         k = 4;
     end
     % atribui a estrtura o valor relacionado a magnetizacao para o H
     struc{1,k} =  f*aa{tt}.ILP;
     % atribui a estrtura o valor relacionado a magnetizacao para o C
     struc{2,k} =  bb{tt}.ILP;
     if k == 4
         % Confiriri os pulsos da sequencia de aquisi??o 
         roh{nto} = tomography({'+II','+YI','+XI','+XX'},struc,2);

         % Deviation matrix
         adesv{nto} = roh{nto};
      
         % Usa a normalizacao do |0,0>
         % Estado normalizado
         roh{nto} = (normalizacao2q).*roh{nto} + eye(4)/4;

         nto = nto + 1;
     end
end

rho_out = roh; 
desv_out = adesv; 

end


%
function [rho_out, desv_out] = tomografia_1q_mult(aa,ff1,t2,norma1q,flag)
% fit da parte A
for k=1:length(aa)
    aa{k} =  fitespec(aa{k},ff1,t2,[0.05 0.2],[-200 200]);
    if flag == 1 
        plotfit(k+10,aa{k},[-200 200]); % para ver todos os espectros
    end
end

% Inicializa contador para fazer mais uma tomografia
nto = 1; % 
f = 1;

for tt = 1:1:length(aa)     
     k = mod(tt,2);
     if k == 0
         k = 2;
     end
     % atribui a estrtura o valor relacionado a magnetizacao para o H
     struc{1,k} =  f*aa{tt}.ILP;
%      % atribui a estrtura o valor relacionado a magnetizacao para o C
%      struc{2,k} =  bb{tt}.ILP;
     if k == 2
         % Confiriri os pulsos da sequencia de aquisi??o 
         roh{nto} = tomography({'+II','+YI'},struc,1);
         
         normalizacao1q = norma1q;
         % Deviation matrix
         adesv{nto} = roh{nto};
      
         % Estado normalizado
         roh{nto} = (normalizacao1q).*roh{nto} + eye(2)/2;

         nto = nto + 1;
     end
end

% % Alternativamente
% % Magnetizacao do Hidrogenio 
% mxh = abs(aa{k-1}.ILP(1) + aa{k-1}.ILP(3))/2; % II
% myh = abs(aa{k-1}.ILP(2) + aa{k-1}.ILP(4))/2; % II 
% mzh = abs(aa{k}.ILP(1) + aa{k}.ILP(3))/2;     % YI         
% rohalt{nto} = (mxh*sigmax + myh*sigmay + mzh*sigmaz)

rho_out = roh; 
desv_out = adesv; 

end

