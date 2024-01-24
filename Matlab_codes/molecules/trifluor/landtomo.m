clc
clear all;

spectropar;
trifluor;


pulso1{1} = decompulse({{'isech200.RF'}},'F1180A',0.3e-3,1,180,{{'n'}});
pulso1{2} = decompulse({{'isech200.RF'}},'F2180A',0.3e-3,2,180,{{'n'}});
pulso1{3} = decompulse({{'isech200.RF'}},'F3180A',0.2e-3,3,180,{{'n'}});
pulso1{4} = decompulse({{'isech200.RF'}},'F190A',0.3e-3,1,90,{{'n'}});
pulso1{5} = decompulse({{'isech200.RF'}},'F290A',0.3e-3,2,90,{{'n'}});
pulso1{6} = decompulse({{'isech200.RF'}},'F390A',0.2e-3,3,90,{{'n'}}); 
pulso1{7} = decompulse({{'isech200.RF'}},'F3aA',0.2e-3,3,140.0,{{'n'}});
pulso1{8} = decompulse({{'isech200.RF','isech200.RF'}},'F1390bA',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso1{9} = decompulse({{'isech200.RF','isech200.RF'}},'F1390cA',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso1{10} = decompulse({{'isech200.RF','isech200.RF'}},'F1390dA',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso1{11} = decompulse({{'isech200.RF','isech200.RF'}},'F1390aA',0.2e-3,[1 3],[90 90],{{'n','n'}}); 


pulso2{1} = decompulse({{'isech200.RF'}},'F1180B',0.3e-3,1,180,{{'n'}});
pulso2{2} = decompulse({{'isech200.RF'}},'F2180B',0.3e-3,2,180,{{'n'}});
pulso2{3} = decompulse({{'isech200.RF'}},'F3180B',0.2e-3,3,180,{{'n'}});
pulso2{4} = decompulse({{'isech200.RF'}},'F190B',0.3e-3,1,90,{{'n'}});
pulso2{5} = decompulse({{'isech200.RF'}},'F290B',0.3e-3,2,90,{{'n'}});
pulso2{6} = decompulse({{'isech200.RF'}},'F390B',0.2e-3,3,90,{{'n'}}); 
pulso2{7} = decompulse({{'isech200.RF'}},'F3aB',0.2e-3,3,140.0,{{'n'}});
pulso2{8} = decompulse({{'isech200.RF','isech200.RF'}},'F1390bB',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso2{9} = decompulse({{'isech200.RF','isech200.RF'}},'F1390cB',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso2{10} = decompulse({{'isech200.RF','isech200.RF'}},'F1390dB',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso2{11} = decompulse({{'isech200.RF','isech200.RF'}},'F1390aB',0.2e-3,[1 3],[90 90],{{'n','n'}}); 


pulso3{1} = decompulse({{'isech200.RF'}},'F1180C',0.3e-3,1,180,{{'n'}});
pulso3{2} = decompulse({{'isech200.RF'}},'F2180C',0.3e-3,2,180,{{'n'}});
pulso3{3} = decompulse({{'isech200.RF'}},'F3180C',0.2e-3,3,180,{{'n'}});
pulso3{4} = decompulse({{'isech200.RF'}},'F190C',0.3e-3,1,90,{{'n'}});
pulso3{5} = decompulse({{'isech200.RF'}},'F290C',0.3e-3,2,90,{{'n'}});
pulso3{6} = decompulse({{'isech200.RF'}},'F390C',0.2e-3,3,90,{{'n'}}); 
pulso3{7} = decompulse({{'isech200.RF'}},'F3aC',0.2e-3,3,140.0,{{'n'}});
pulso3{8} = decompulse({{'isech200.RF','isech200.RF'}},'F1390bC',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso3{9} = decompulse({{'isech200.RF','isech200.RF'}},'F1390cC',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso3{10} = decompulse({{'isech200.RF','isech200.RF'}},'F1390dC',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso3{11} = decompulse({{'isech200.RF','isech200.RF'}},'F1390aC',0.2e-3,[1 3],[90 90],{{'n','n'}}); 


pulso4{1} = decompulse({{'isech200.RF'}},'F1180D',0.3e-3,1,180,{{'n'}});
pulso4{2} = decompulse({{'isech200.RF'}},'F2180D',0.3e-3,2,180,{{'n'}});
pulso4{3} = decompulse({{'isech200.RF'}},'F3180D',0.2e-3,3,180,{{'n'}});
pulso4{4} = decompulse({{'isech200.RF'}},'F190D',0.3e-3,1,90,{{'n'}});
pulso4{5} = decompulse({{'isech200.RF'}},'F290D',0.3e-3,2,90,{{'n'}});
pulso4{6} = decompulse({{'isech200.RF'}},'F390D',0.2e-3,3,90,{{'n'}}); 
pulso4{7} = decompulse({{'isech200.RF'}},'F3aD',0.2e-3,3,140.0,{{'n'}});
pulso4{8} = decompulse({{'isech200.RF','isech200.RF'}},'F1390bD',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso4{9} = decompulse({{'isech200.RF','isech200.RF'}},'F1390cD',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso4{10} = decompulse({{'isech200.RF','isech200.RF'}},'F1390dD',0.2e-3,[1 3],[90 90],{{'n','n'}}); 
pulso4{11} = decompulse({{'isech200.RF','isech200.RF'}},'F1390aD',0.2e-3,[1 3],[90 90],{{'n','n'}}); 


ss1 =  compseq('landauer1At',pulso1);
roh = mol.Iz{1} + mol.Iz{2} + mol.Iz{3}  ;
[roh2 spec]= simseq(ss1,roh,[-200 200]);
seq2spec('land1At',ss1,pulso1);

ss2 =  compseq('landauer1Bt',pulso2);
roh = mol.Iz{1} + mol.Iz{2} + mol.Iz{3}  ;
[roh2 spec]= simseq(ss2,roh,[-200 200]);
seq2spec('land1Bt',ss2,pulso2);

ss3 =  compseq('landauer1Ct',pulso3);
roh = mol.Iz{1} + mol.Iz{2} + mol.Iz{3}  ;
[roh2 spec]= simseq(ss3,roh,[-200 200]);
seq2spec('land1Ct',ss3,pulso3);

ss4 =  compseq('landauer1Dt',pulso4);
roh = mol.Iz{1} + mol.Iz{2} + mol.Iz{3}  ;
[roh2 spec]= simseq(ss4,roh,[-200 200]);
seq2spec('land1Dt',ss4,pulso4);


clc
pp1(1)=1e-3*ss1{end-1}.desc.power*4095/21.097
pp1(2)=1e-3*ss1{end-5}.desc.power*4095/21.097
pp1(3)=1e-3*ss1{end-9}.desc.power*4095/21.097
pp1(4)=1e-3*ss1{end-13}.desc.power*4095/21.097

pp2(1)=1e-3*ss2{end-2}.desc.power*4095/21.097
pp2(2)=1e-3*ss2{end-6}.desc.power*4095/21.097
pp2(3)=1e-3*ss2{end-12}.desc.power*4095/21.097
pp2(4)=1e-3*ss2{end-14}.desc.power*4095/21.097


pp3(1)=1e-3*ss3{end-2}.desc.power*4095/21.097
pp3(2)=1e-3*ss3{end-6}.desc.power*4095/21.097
pp3(3)=1e-3*ss3{end-12}.desc.power*4095/21.097
pp3(4)=1e-3*ss3{end-14}.desc.power*4095/21.097

pp4(1)=1e-3*ss4{end-4}.desc.power*4095/21.097
pp4(2)=1e-3*ss4{end-8}.desc.power*4095/21.097
pp4(3)=1e-3*ss4{end-12}.desc.power*4095/21.097



