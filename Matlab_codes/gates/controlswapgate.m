function U = controlswapgate(nq, controle, alvo1, alvo2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = controlswapgate(nq, controle, alvo1, alvo2)
%nq = number of qubits
%alvo1 = target qubit 1 
%alvo2 = target qubit 2 
%controle = control qubits
%U = controlgate(5,[1 3],2,4) -> U is a 5 qubit unitary that swap qubits 2
%and 4 if  1 and 3 are on state 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    
U = controlgate(nq, [controle alvo1], alvo2, xgate)*...
        controlgate(nq, [controle alvo2], alvo1, xgate)*...
        controlgate(nq, [controle alvo1], alvo2, xgate);
    
