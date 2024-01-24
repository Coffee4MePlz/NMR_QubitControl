function U = controlgate(nq,controle,alvo,Uop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = controlgate(nq,controle,alvo,Uop);
%nq = number of qubits
%alvo = target qubits 
%controle = control qubits
%Uop = one qubit gate
%U = controlgate(5,[1 3],[2 4],Uop) -> U is a 5 qubit unitary of a control gate
%where Uop is applied to qubits 2 and 4  if 1 and 3 are on state 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

alvo = sort(alvo);

controle = sort([controle]);

na = length(alvo);

U = 1;

Up = [1 0; 0 0];

Dw = [0 0;0 1];

aux(1:nq) = 1;    aux(controle) = 0;

for m=1:na
    Uca = Uop;
    
    for k=alvo(m)-1:-1:1;     
        Uca = kron(Dw,Uca) + kron(Up,Uca^aux(k)); 
    end 

    for k=alvo(m)+1:1:nq;  
        Uca = kron(Uca,Dw) + kron(Uca^aux(k),Up); 
    end 
    
    U = Uca*U;

end

end

