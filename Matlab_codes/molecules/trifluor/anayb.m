clc

clear all

spectropar;

trifluor;

[freq1 tran ff1 opR opI] = transitions(mol,2);
  
ref= varian2mat('swap',64*1024);
 
ref= apodize(ref);
  
ref{1} = fitespec(ref{1},ff1,mol.T2(2)/pi,[0.1 1],[-400 400]);
  
 plotfit(1,ref{1},[-200 200]); pause(0.2);
  
 mref = sum(ref{1}.ILP(3));
  
  
 %ss= varian2mat('oneyb5sat',64*1024);
 %ss= varian2mat('oneyb6notsat',64*1024);
%  ss= varian2mat('oneybnotsattwox',64*1024);
ss= varian2mat('oneyb7A2Bnotsat',64*1024);
 
 ss= apodize(ss);
 
 
 for k=1:length(ss);
 
  ss{k} = fitespec(ss{k},ff1,mol.T2(2)/pi,[0.1 1],[-400 400]);
  
   plotfit(1,ss{k},[-200 200]); pause(0.2);
  
   mm7(k) = sum(ss{k}.ILP(1:2:8));
 end
 
 mm7/mref
 
 
 
 plot(mm7/mref)
% 
% plot(mm3/mref)
% 
% plot(mm4/mref)










