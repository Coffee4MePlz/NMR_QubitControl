
close all
clear all
% loc = pwd;
% cd /home/amsouza/matlab/;
% addpath([pwd '/others/']) 
% addpath([pwd '/nmr/']) 
% addpath([pwd '/gates/']) 
% addpath([pwd '/shapes/']) 
% addpath([pwd '/exps/'])
% addpath([pwd '/control/'])
% eval(['cd ' loc]);

clc

spectropar;
trifluor;


grape.nome = 'F290';
grape.U = multigate(3,2,rotx(pi/2));
grape.np = 400;
grape.temp = 1000e-6;
grape.maxpower = 5e3;
grape.chanel = [1];
grape.rf = [0];
grape.chem{1} = [];
grape.chem{2} = [];
grape.chem{3} = [];

findgrape(mol,spectro,grape);
unix(['cp ' grape.nome 'Obs.RF' ' /home/amsouza/matlab/shapes/' grape.nome 'Obs.RF'])


pulso{1} = decompulse({{'F290Obs.RF'}},'F290',1e-3,2,90,{{'g'}});

 ss =  compseq('exemplo2test',pulso,'nome');
 roh = mol.Iz{1} + mol.Iz{2} + mol.Iz{3}  ;
 [roh2 spec]= simseq(ss,roh,[-200 200],2);