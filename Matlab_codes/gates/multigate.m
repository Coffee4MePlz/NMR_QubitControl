function U = multigate(nq,alvo,Uop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = controlgate(nq,alvo,Uop);
%nq = number of qubits
%alvo = target qubits 
%Uop = one qubit gate
%U = controlgate(5,[2 4],Uop) -> U is a 5 qubit unitary where Uop is applied to 
%qubits 2 and 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

aux = cell(1,nq);

for jj=1:nq; aux{jj} = eye(2); end

for jj=1:length(alvo); aux{alvo(jj)} = Uop; end
   
U  = aux{1}; 

for k=1:(nq -1); U = kron(U,aux{k+1}); end

end
