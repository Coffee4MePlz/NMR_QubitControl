function f=errorbarmat(r1,error) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Roberto Serra
% f=errorbarmat(r1,error)
%
% similate random error in an n qubit density matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nqubit = length(r1);
if nqubit == 2
    a = rand;
    b = rand;
    c = rand + 1i*rand;
    roh1 = [ a c; conj(c) b];
    rohq = roh1/real(trace(roh1));
else
    nqubit = log2(nqubit);
    for n = 1:nqubit
        a = rand;
        b = rand;
        c = rand + 1i*rand;
        roh1 = [ a c; conj(c) b];  % hermitian random 2x2 matrix
        roh1 = roh1/real(trace(roh1)); % normalized hermitian random 2x2 matrix
        roh{n} = roh1;
        if n == 1
           rohb =  roh{n};
        else    
            roha = roh{n};
            rohq = kron(roha,rohb); % nxn normalized hermitian random 2x2 matrix
            rohb = rohq;
        end
    end
end  

x = randi([0 1]);
aux = r1 + ((-1)^x)*error*rohq;
aux = aux/abs(trace(aux));% matrix normalization

% Max likellyhood - need to be implemented.
%

f = aux;

    

