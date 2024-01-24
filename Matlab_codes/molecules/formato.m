
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definicoes Molecula %%%%%%%%%%%%%%%%%%%%%%%%%%
 
global mol;
 
%NOME DA MOLECULA
mol.nome = 'formato';
 
%NUMERO DE SPINS
mol.nspin = 2;
 
%RANK DE CADA SPIN
mol.stipo = [1/2 1/2];
 
%T2 
mol.T2 = [0.232; 0.317; 0.05];
 
%DESLOCAMENTO QUIMICO DE CADA SPIN % na Ressonacia ? zero
mol.dq = [0  0];
 
%ACOPLAMENTO J (escalar) - atualiazar a cada exp.
mol.j = [0    194.65;
         0     0];
 
%Tipo do acoplamento (0 -> ZZ e 1 -> XX + YY + ZZ)
mol.jtipo = [0 0;
             0 0];
 
%Operadores de Spin - Em unidades de h/ 2*pi = ??
%Forma de acesso pela comand window: mol.Ix{1} e etc.
%Ex.: mol.Ix{2} = kron(I,Sx)
[mol.Ix mol.Iy mol.Iz mol.Ich] = GenPauli(mol.stipo);
 
%Hamiltoniano - GenH
[mol.Hzee mol.Hint] = GenH(mol);
