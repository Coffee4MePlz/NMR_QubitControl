function [spec]  = specgen(roh,chanel)
% Generate spectrum
%[spec]  = specgen(roh,chanel)
% roh -> input density matrix
% chanel -> chanel to caculate the spectrum
% spec - matlab structure with the simulate FID and spectrum 

global mol spectro

%Encontra autovalores e autovetores da Hamiltoiniana
[eigvec eigval] = eig(mol.Hint + mol.Hzee);

eigval = diag(eigval);

%Matriz densidade e operadores na base H
roh = (eigvec')*roh*(eigvec);

for n=1:length(spectro.chanel{chanel})
Ixy{n} = (eigvec')*(mol.Ix{spectro.chanel{chanel}(n)}+ 1i*mol.Iy{spectro.chanel{chanel}(n)})*(eigvec); 
end

%Tamanho do fid
nfid  = spectro.nfid;
sw = spectro.sw;
spec.fid = zeros(1,nfid);
t = (0:nfid-1)*(1/sw);


%Calcula fid espectro
for k=1:prod(2*mol.stipo +1);

    for m=1:prod(2*mol.stipo +1);

          fr = (eigval(m) - eigval(k));
    
          if fr < 2*pi*sw/2 & fr > - 2*pi*sw/2
          
                for n=1:length(spectro.chanel{chanel})
     
                    spec.fid = spec.fid + roh(k,m)*Ixy{n}(m,k)*...
                    exp(1i*fr*t -t/(mol.T2(spectro.chanel{chanel}(n))));
                end
          end
    
    end 
   
end

spec.esp = fftshift(fft(spec.fid,spectro.nesp));

spec.sw = sw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

