function matrix = ztrace(roh, n);
%
%
%    ztrace(roh, qubit) = traço parcial de roh relativo ao qubit 
%
%    qubit = 1, 2,..., n; onde roh tem dimensão 2^n x 2^n;  
%
% 

[N, N] = size(roh);

numero = N / (2^n);  qmax = log2(N);

vet = 1:1:N;


ind1 = []; ind2 = [];


for jj = 1:1:(N / 2^(qmax - n + 1));

   
      ind1 = [ind1 (1 + (jj - 1) * 2 * numero):1:(numero + (jj - 1) * 2 * numero)]; 
      
      ind2 = [ind2 (1 + (jj - 1) * 2 * numero + numero):1:(2 * numero + (jj - 1) * 2 * numero)];   
   
end;
   
matrix = roh(ind1, ind1) + roh(ind2, ind2);
   
