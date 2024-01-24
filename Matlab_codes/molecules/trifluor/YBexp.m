clc
clear all;

%spectrometer configuration
spectropar;

%define the molecule
trifluor;

%decompose the pulse used in the sequence. All information abiut the pulses
%is wirtten in the variable pulse
pulso{1} = decompulse({{'isech200.RF'}},'F1180',0.3e-3,1,180,{{'n'}});
pulso{2} = decompulse({{'isech200.RF'}},'F2180',0.3e-3,2,180,{{'n'}});
pulso{3} = decompulse({{'isech200.RF'}},'F3180',0.2e-3,3,180,{{'n'}});
pulso{4} = decompulse({{'isech200.RF'}},'F190',0.3e-3,1,90,{{'n'}});
pulso{5} = decompulse({{'isech200.RF'}},'F290',0.3e-3,2,90,{{'n'}});
pulso{6} = decompulse({{'isech200.RF'}},'F390',0.2e-3,3,90,{{'n'}}); 
pulso{7} = decompulse({{'isech200.RF','isech200.RF'}},'F2390',0.2e-3,[2 3],[90 90],{{'n','n'}});

%compilation of exemplo1test
 ss =  compseq2('YB',pulso,'yb');

 %sumilation  
 roh = mol.Iz{1} + mol.Iz{2} + mol.Iz{3}  ;
 [roh2 spec]= simseq(ss,roh,[-200 200],1);

 %the compilded sequence and pulses are passed to the spectrometer 
 %seq2spec('seqtest',ss,pulso);