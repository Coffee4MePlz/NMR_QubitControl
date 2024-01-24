
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definicoes Molecula %%%%%%%%%%%%%%%%%%%%%%%%%%

global mol;

%NOME DA MOLECULA
mol.nome = 'trifluor';

%NUMERO DE SPINS
mol.nspin = 3;

%RANK DE CADA SPIN
mol.stipo = [1/2 1/2 1/2];

%T2
mol.T2 = [0.2166; 0.3436; 0.118];

%DESLOCAMENTO QUIMICO DE CADA SPIN %Hz
mol.dq = [11833.8 0 -17324.1] - 11833.8;

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

