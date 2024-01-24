
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definicoes Molecula %%%%%%%%%%%%%%%%%%%%%%%%%%

global mol;
 
%NOME DA MOLECULA
mol.nome = 'glicina';
 
%NUMERO DE SPINS
mol.nspin = 4;
 
%RANK DE CADA SPIN
mol.stipo = [1/2 1/2 1/2 1/2];
 
mol.T2 = [0.1518; 0.1532; 0.47; 0.207];

%DESLOCAMENTO QUIMICO DE CADA SPIN
mol.dq = [0 0 [[9748.7 -6715.7]+ 6715.7]];

%ACOPLAMENTO J 
mol.j = [0    0   5.3   143.6;
         0    0   5.3   143.6;
         0    0   0     53.6;
         0    0   0      0];
     
%Tipo do acoplamento (0 -> ZZ e 1 -> XX + YY + ZZ)
mol.jtipo = [0 0 0 0;
             0 0 0 0;
             0 0 0 0;
             0 0 0 0];
 

%Operadores de Spin


 [mol.Ix mol.Iy mol.Iz mol.Ich] = GenPauli(mol.stipo);

 [mol.Hzee mol.Hint] = GenH(mol);

