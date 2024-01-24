function [Ix Iy Iz Ich] = GenPauli(Stipo)
%Operadores de Spin 

%global - Declare variables as global
global spectro;

% number of qubits
nq = length(Stipo);

ss = 2*Stipo+1;

%cell(x,y) - create a cell array with dimensions "x by y"
Ix = cell(1,nq);
Iy = cell(1,nq);
Iz = cell(1,nq);

for k=1:nq
    
    %Programa "mat_ixyz" é responsável pelo cálculo dos operadores de spin,
    %com base no valor de spin inserido em "cloroformio".
    %O resultado é uma "cell array", os elementos são os respectivos operadores.
    [sx sy sz] = mat_ixyz(Stipo(k));
   
    %Note: In a script file which contains commands and function definitions. 
    %Functions must be at the end of the file. Script files cannot have the same name 
    %as a function in the file.
   
    %Calculo dos produtos tensoriais -  Função local "mykron".
    Ix{k} = (mykron(nq,k,sx,ss));
    Iy{k} = (mykron(nq,k,sy,ss));
    Iz{k} = (mykron(nq,k,sz,ss));
end   

Ich = cell(2,2);
%Ich será uma matriz que abriga em seus elementos, alguns produtos tensoriais
%calculados anteriormente.
%Por exemplo: mol.Ich{1,1} = mol.Ix{1} = kron(Sx,I)
for k=1:2
    Ich{k,1} = zeros(prod(ss));
    Ich{k,2} = zeros(prod(ss));
    for m=1:length(spectro.chanel{k});
        Ich{k,1} = Ich{k,1} + Ix{spectro.chanel{k}(m)};
        Ich{k,2} = Ich{k,2} + Iy{spectro.chanel{k}(m)};
    end
    
end

end


%Essa função é responsável por criar os produtos tensoriais do tipo
%kron(I,X) ou kron(X,I)
%Forma de acesso pela comand window: mol.Ix{1} e etc.
%mol.Ix{2} = kron(I,Sx) and  %mol.Ix{1} = kron(Sx,I)
%mol.Iy{2} = kron(I,Sy) and  %mol.Iy{1} = kron(Sy,I)
%mol.Iz{2} = kron(I,Sz) and  %mol.Iz{1} = kron(Sz,I)
function gat = mykron(nq,na,Unit,ss);
aux = cell(1,nq);

for jj=1:nq; aux{jj} = eye(ss(jj)); end

for jj=1:length(na); aux{na(jj)} = Unit; end
   
gat  = aux{1}; for k=1:(nq -1); gat = kron(gat,aux{k+1}); end

end