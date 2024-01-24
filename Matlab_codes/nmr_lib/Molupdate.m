function [tof] = Molupdate(data,spin,mol,guesstal,guessf);

ns = length(spin);

dq = mol.dq;

mol2 =mol;

for k=1:ns
	   mol2.dq = mol.dq - mol.dq(spin(k)); 
	   
	   [freq1 tran ff1 opR opI] = transitions(mol2,spin(k));
       
	   ss= varian2mat(data{k},64*1024);
       
	   ss= apodize(ss);
       
	   ss{1} = fitespec(ss{1},ff1+ guessf(k),guesstal(k),[guesstal(k) 5],[-400 400]);
       
	   plotfit(k,ss{1},[-200 200]); pause(0.2);	
	   
	   df(k) = ss{1}.INLP(2);
	   
	   dq = mol.dq(spin(k)) + df(k) + guessf(k);
	   
	   tof(k) = ss{1}.tof + df(k)+ guessf(k);
	   
	   tal(k) = ss{1}.INLP(1);
end

disp(df)

disp(tof)

disp(tal)

