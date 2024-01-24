function U = swapgate(nq,alvo1,alvo2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = swapgate(nq,alvo1,alvo2);
%nq = number of qubits
%alvo1 = target qubit 1 
%alvo2 = target qubit 2 
%
%U = controlgate(5,2,4) -> U is a 5 qubit unitary that swap qubits 2
%and 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


U = controlgate(nq,alvo1,alvo2,xgate)*...
     controlgate(nq,alvo2,alvo1,xgate)*...
     controlgate(nq,alvo1,alvo2,xgate);
