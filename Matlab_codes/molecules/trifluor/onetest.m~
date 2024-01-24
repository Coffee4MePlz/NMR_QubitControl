clc

clear all;

spectropar;

trifluor;

pulso{1} = decompulse({{'isech200.RF'}},'F1180',0.3e-3,1,180,{{'n'}});
pulso{2} = decompulse({{'isech200.RF'}},'F2180',0.25e-3,2,180,{{'n'}});
pulso{3} = decompulse({{'isech200.RF'}},'F3180',0.2e-3,3,180,{{'n'}});
pulso{4} = decompulse({{'isech200.RF'}},'F190',0.32e-3,1,90,{{'n'}});
pulso{5} = decompulse({{'isech200.RF'}},'F290',0.24e-3,2,90,{{'n'}});
pulso{6} = decompulse({{'isech200.RF'}},'F390',0.21e-3,3,90,{{'n'}}); 

pulso{7} = decompulse({{'F190AYBAObs.RF'}},'F190g',1.0e-3,1,90,{{'g'}}); 
pulso{8} = decompulse({{'F290AYBAObs.RF'}},'F290g',1.0e-3,2,90,{{'g'}}); 
pulso{9} = decompulse({{'F390AYBAObs.RF'}},'F390g',1.0e-3,3,90,{{'g'}}); 


 ss =  compseq2('oneYB',pulso,'oneYB7');
 
 
 roh = mol.Iz{1} + mol.Iz{2} + mol.Iz{3}  ;
 [roh2 spec]= simseq(ss,roh,[-200 200],1);
 seq2spec(['oneYB7'],ss,pulso);

  
  
  
  
 