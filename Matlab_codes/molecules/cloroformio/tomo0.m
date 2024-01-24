
clear all

spectropar; % Leitura de parametros  no arquivo arquivo spectropar.m

cloroformio; % Leitura dos parametros da molecula

%Find the frequncies and transtions of the spin "sp" of the molecule "mol"
[tran1 tran1 ff1 opr] = transitions(mol,1);

[tran2 tran2 ff2 opr] = transitions(mol,2);

 % Referencia para o Hidrogenio
 ff1 = ff1+0.1; %atualizar para ver se esta fora da ressoancia 

  % atualizar'referenciah' para o nome do arquivo do fid
  % Importa os dados do experimento
  %spechr= varian2mat('referenciah',16*1024); 
  spechr= varian2mat('refh90',128*1024);  
  %spechr= varian2mat('calib180h',32*1024);
  %spechr= varian2mat('refh2',32*1024);
  spechr{1} = fitespec(spechr{1},ff1,0.67,[0.01 0],[-200 200]);
  plotfit(1,spechr{1},[-200 200]); 
  
  % Referencia para o Carbono
  
  %ff2 = ff2-0.4;%atualizar
  
%   speccr= varian2mat('referenciac',16*1024);
%    speccr= varian2mat('refc',64*1024);
%    speccr{1} = fitespec(speccr{1},ff2,0.175,[0.001 0.0],[-200 200]);
%    plotfit(1,speccr{1},[-200 200]); 
%  
%  nh = abs(spechr{1}.ILP(1) + spechr{1}.ILP(3))/2;
%  nc = abs(speccr{1}.ILP(1) + speccr{1}.ILP(3))/2;
%  f =4*nc/nh;
 
 
 % colocar o nome do arquivo sem o .ifd
 % bb=varian2mat('calb180h',64*1024);
 % bb=varian2mat('pseunovo',64*1024);
 % bb=varian2mat('terh',64*1024);
   bb=varian2mat('ib1fz1z2',128*1024);
   
   
 %atualizar tempo

 for k=1:length(bb); 	
 	bb{k} =  fitespec(bb{k},ff1,0.67,[0.6 0.6],[-200 200]);
     plotfit(k,bb{k},[-200 200]); 
 end
 

 %est1 = [0.783 0; 0 0.217]; 

	
 est1 = [1 0; 0 0]; 

 
 nto = 1;
 for tt = 1:1:length(bb)
     
     k = mod(tt,2);
     
     if k == 0
         k = 2;
	 end
     
	 
     	 

     struc{1,k} =  bb{tt}.ILP;

     if k == 2
         
         rohc{nto} = tomography({'+II','+XI'},struc,2,[1]);
         
         %normalizacao1q = ((1.0/2.0)/rohc{1,1}(1,1));
         
        normalizacao1q = 5.016004197205741e-06
         rohc{nto} = (normalizacao1q).*rohc{nto} + eye(2)/2.
         
         fidelmat(rohc{nto},est1)

		 showmat(5,rohc{nto}); pause(0.2);
         
         

         nto = nto+1
     end
 
     % Calcula a magnetizao em X a partir do expectro fitado
     mx(tt) = (bb{tt}.ILP(1) + bb{tt}.ILP(3))/2.;

     % Calcula a magnetizao em Y a partir do expectro fitado
     my(tt) = (bb{tt}.ILP(2) + bb{tt}.ILP(4))/2.;
     
 eixox(tt) =  tt
 end
 %figure(100)
 %plot(eixox,mx);
 %figure(200)
 %plot(eixox,my);
 
% for k=1:38
%    plotfit(k,specx{1}{k},[-200 200]); 
% end

  save('tomografia2parte3novo.mat','rohc');
