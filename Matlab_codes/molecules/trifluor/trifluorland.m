
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definicoes Molecula %%%%%%%%%%%%%%%%%%%%%%%%%%

global mol;

%NOME DA MOLECULA
mol.nome = 'trifluor';

%NUMERO DE SPINS
mol.nspin = 3;

%RANK DE CADA SPIN
mol.stipo = [1/2 1/2 1/2];

%T2
mol.T2 = [0.0757    0.0897    0.0760];

%DESLOCAMENTO QUIMICO DE CADA SPIN
mol.dq = [9839.6 -1995.3 -19320.1] +1995.3 ;
% mol.dq = [0  -11835.8  -29162.7];
%ACOPLAMENTO J 
mol.j = [0   69.86  47.65;
         0    0   -128.1;
         0    0   0];

%Tipo do acoplamento (0 -> ZZ e 1 -> XX + YY + ZZ)
mol.jtipo = [0 0 0;
             0 0 0;
             0 0 0];

%Operadores de Spin


 [mol.Ix mol.Iy mol.Iz mol.Ich] = GenPauli(mol.stipo);

 [mol.Hzee mol.Hint] = GenH(mol);

