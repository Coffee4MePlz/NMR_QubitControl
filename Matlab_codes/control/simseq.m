
function [roh  spec]= simseq(struc,roh, opt, nfig);
% Simulate the compiled sequence "struc"
% on the iniital density matrix roh


global mol spectro

roh2 = roh;  m =1;
roh2= struc{end}.Uid{m}*roh*struc{end}.Uid{m}';

U = 1;
 for k=1:length(struc)-1

    switch struc{k}.tip
               case 'pulse'

                      nchan  = length(struc{k}.desc.chan);
                      
                      if nchan == 1;
                            if struc{k}.desc.chan{1} == 1;
                          [Ur roh]=shapedpulse(struc{k}.desc.formaspec{1},1e6*struc{k}.desc.pw,0,0,struc{k}.desc.power,roh);
                            else 
                             [Ur roh]=decshapedpulse(struc{k}.desc.formaspec{1},1e6*struc{k}.desc.pw,0,0,struc{k}.desc.power,roh);
                            end
                      else
                              if struc{k}.desc.chan{1} == 1; n1 =1; n2 = 2; else n1 = 2; n2 =1; end;
                             [Ur roh]=simshapedpulse(struc{k}.desc.formaspec{n1},struc{k}.desc.formaspec{n2},...
                             1e6*struc{k}.desc.pw,0,0,0,0,struc{k}.desc.power(n1),struc{k}.desc.power(n2),roh);
                      end
                      
               case 'delay'
               dt = struc{k}.desc; [U roh]=delay(dt,roh);
              
               case 'grad'
               roh = diag(diag(roh));
			   roh2 = diag(diag(roh2));
               m=m+1; roh2= struc{end}.Uid{m}*roh2*struc{end}.Uid{m}';
               
      end
 end


      for m=1:mol.nspin;
      Uc = multigate(mol.nspin,m,rotz(pi*struc{end}.rz(m)/180));
      roh = Uc*roh*Uc';
end

[spec]  = specgen(roh,1);


if nargin == 2; plotesp(2,spec); else plotesp(nfig,spec,opt); end

for m=1:mol.nspin;
     Uc = multigate(mol.nspin,m,rotphi(pi*struc{end}.virtua(m),struc{end}.phasevirtua(m)*pi/180));
      roh = Uc*roh*Uc';
end

showmat(10,roh)
disp(fidelmat(roh2,roh))




