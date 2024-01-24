
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definicoes Molecula %%%%%%%%%%%%%%%%%%%%%%%%%%

global mol;

global mol;
 
%NOME DA MOLECULA
mol.nome = 'tricloro';
 
%NUMERO DE SPINS
mol.nspin = 3;
 
%RANK DE CADA SPIN
mol.stipo = [1/2 1/2 1/2];
 
%T2
mol.T2 = [0.2; 0.1; 0.05];
 
%DESLOCAMENTO QUIMICO DE CADA SPIN
mol.dq = [0  -415 415];
 
%ACOPLAMENTO J 
mol.j = [0   201 8;
         0    0  103;
         0    0   0];
 
%Tipo do acoplamento (0 -> ZZ e 1 -> XX + YY + ZZ)
mol.jtipo = [0 0 0;
             0 0 1;
             0 0 0];
 

%Operadores de Spin


 [mol.Ix mol.Iy mol.Iz mol.Ich] = GenPauli(mol.stipo);

 [mol.Hzee mol.Hint] = GenH(mol);

