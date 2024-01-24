clear all

spectropar; % Leitura de parametros  no arquivo arquivo spectropar.m

cloroformio; % Leitura dos parametros da molecula

%Find the frequncies and transtions of the spin "sp" of the molecule "mol"
[tran1 tran1 ff1 opr] = transitions(mol,1);

[tran2 tran2 ff2 opr] = transitions(mol,2);

  
  % Referencia para o Carbono
  
   ff2 = ff2 ;%atualizar
  
%   speccr= varian2mat('referenciac',16*1024);
    speccr= varian2mat('refc',64*1024);
    speccr{1} = fitespec(speccr{1},ff2,0.175,[0.001 0.0],[-200 200]);
    plotfit(1,speccr{1},[-200 200]); 