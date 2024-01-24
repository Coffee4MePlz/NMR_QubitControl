function U = zzgate(nq,ang,alvo1,alvo2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = zzgate(nq,ang,alvo1,alvo2);
%nq = number of qubits
%alvo1 = target qubit 1 
%alvo2 = target qubit 2 
%ang = angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ZZ = multigate(nq,[alvo1 alvo2],zgate)/4;

U = expm(-1i*ang*ZZ);
