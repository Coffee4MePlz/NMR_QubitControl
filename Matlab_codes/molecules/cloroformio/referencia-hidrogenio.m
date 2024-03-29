clear all

spectropar; % Leitura de parametros  no arquivo arquivo spectropar.m

cloroformio; % Leitura dos parametros da molecula

%Find the frequncies and transtions of the spin "sp" of the molecule "mol"
[tran1 tran1 ff1 opr] = transitions(mol,1);

[tran2 tran2 ff2 opr] = transitions(mol,2);

 % Referencia para o Hidrogenio
  ff1 = ff1; %atualizar para ver se esta fora da ressoancia 

  % atualizar'referenciah' para o nome do arquivo do fid
  % Importa os dados do experimento
  spechr= varian2mat('refh',64*1024);  
  spechr{1} = fitespec(spechr{1},ff1,0.505,[0.001 0],[-200 200]);
  plotfit(1,spechr{1},[-200 200]); 