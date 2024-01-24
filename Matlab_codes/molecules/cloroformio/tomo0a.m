

clear all

spectropar;

cloroformio;

[tran1 tran1 ff1 opr] = transitions(mol,1);

[tran2 tran2 ff2 opr] = transitions(mol,2);

ff1 = ff1+1.5;%atualizar


  spechr= varian2mat('referenciah',16*1024);
  %spechr= varian2mat('refh2',32*1024);
  spechr{1} = fitespec(spechr{1},ff1,0.33,[0.2 3],[-200 200]);
  plotfit(1,spechr{1},[-200 200]); 
  
  
  ff2 = ff2-0.4;%atualizar
  
%   speccr= varian2mat('referenciac',16*1024);
%   %speccr= varian2mat('refc',32*1024);
%   speccr{1} = fitespec(speccr{1},ff2,0.30,[0.2 3],[-200 200]);
%   plotfit(1,speccr{1},[-200 200]); 
%  
%  nh = abs(spechr{1}.ILP(1) + spechr{1}.ILP(3))/2;
%  nc = abs(speccr{1}.ILP(1) + speccr{1}.ILP(3))/2;
%  f =4*nc/nh;
 
 

  bb=varian2mat('nome',16*1024);


 

 for k=1:length(bb); 	
 	bb{k} =  fitespec(bb{k},ff1,0.30,[0.2 3],[-200 200]);
     plotfit(1,bb{k},[-200 200]); 
 end
 

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
         
         normalizacao1q = ((1.0/2.0)/rohc{1,1}(1,1));

         rohc{nto} = (normalizacao1q).*rohc{nto} + eye(2)/2.
         
%          fidelmat(rohc{nto},est2)

		 %showmat(5,rohc{nto}); pause(0.2);

         nto = nto+1
     end
 
 
 
  end
  save('tomografia2parte3novo.mat','rohc');
