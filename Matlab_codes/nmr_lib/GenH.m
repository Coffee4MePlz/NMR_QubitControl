function  [Hz Hj] = GenH(mol);

Hj = zeros(size(mol.Iz{1}));

%%%%%%%%%%%%%%%%%%%%%% Interaction Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%

%O primeiro laço define se temos apenas acoplamento escalar ou outros tipos
%de acoplamento, entre os spins. (0 -> ZZ e 1 -> XX + YY + ZZ)
  for k=1:mol.nspin; 
        for m=(k+1):mol.nspin;
               
               if mol.jtipo(k,m) == 1;
                 
                Hj = Hj + 2*pi*mol.j(k,m)* ...
               (mol.Ix{k}*mol.Ix{m} + mol.Iy{k}*mol.Iy{m} + mol.Iz{k}*mol.Iz{m});
               end
               
               %2 spins; heteronucleares; C e H.
               %mol.j(k,m) - contante de cacoplamento escalar
               %mol.Iz{k}*mol.Iz{m} - equivale ao produto tensorial sigmaZ
               %SigmaZ.
               if mol.jtipo(k,m) == 0;
               Hj = Hj + 2*pi*mol.j(k,m)*mol.Iz{k}*mol.Iz{m};
               end
        end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%% Zeeman Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hz = 0;

for k=1:mol.nspin;
    Hz = Hz + 2*pi*mol.dq(k)*mol.Iz{k};
end