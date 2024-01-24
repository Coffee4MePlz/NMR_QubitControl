clc
clear all

spectropar;
trifluor;

 [freq1 tran ff1 opR opI] = transitions(mol,2);
  %ref= varian2mat('reference',64*1024);
  %ref= apodize(ref);
  %ref{1} = fitespec(ref{1},ff1,0.98,[3 4],[-400 400]);
  %plotfit(1,ref{1},[-200 200])
  
  %mref = sum(ref{1}.ILP(1:2:8))
  
   ss= varian2mat('exp01',64*1024);
   ss= apodize(ss);
  
  clear mx my
  for k=1:length(ss)
 	 ss{k} = fitespec(ss{k},ff1,0.98,[1 1],[-400 400]);
	 
	 %ss{k} = fitespec(ss{k},ff1,1./(pi*mol.T2(2)),[1./(pi*mol.T2(2)) 3],[-400 400]);
     plotfit(1,ss{k},[-200 200])
     pause(0.1)
     mx(k) = sum(ss{k}.ILP(1:2:8)); %/mref;
     my(k) = sum(ss{k}.ILP(2:2:8)); %/mref;
  end
  
  
  plot(1:length(ss),(mx),'ob',1:length(ss),(my),'or')
 
 