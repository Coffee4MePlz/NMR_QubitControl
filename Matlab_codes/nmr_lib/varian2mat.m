function struc = varian2mat(arq,SI,nor)
%read the varian file
%struc = varian2mat(arq,SI,nor)
% arq - > name of file 
%SI - > number of points to make FFT (fast fourier transform)
%nor -> normalization factor (optional)
%struc - > matlab structure with data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROCPAR
fclose all;
% Para windows atualizar conforme a pasta dos fids e o sistema operacional
% -------->    atualizar tambem o caminho na linha 35   <-----------
%id=fopen([ 'C:\Users\John\Documents\MATLAB\dadosfid\' arq '.fid\procpar'], 'r', 'b');
% Para Linux espectometro
id=fopen([ '/home/user/nmr/fids/' arq '.fid/procpar'], 'r', 'b');


aux2 = [];

cont  = 1;
while cont == 1
aux1 = fgetl(id);
if aux1 == -1; break; end
aux2 = str2mat(aux2,deblank(aux1));
end

eval(['SW =' aux2(strmatch('sw'  , aux2)+1,2:end) ';']);
eval(['rp =' aux2(strmatch('rp '  , aux2)+1,2:end) ';']);
eval(['tof =' aux2(strmatch('tof '  , aux2)+1,2:end) ';']);
eval(['lb =' aux2(strmatch('lb '  , aux2)+1,2:end) ';']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FID

fclose all;
id = fopen([ '/home/user/nmr/fids/' arq '.fid/fid'],'r','b');			
[a, count] = fread(id, 6, 'int32');
nblocks = a(1);
ntraces = a(2);
np = a(3);
ebytes = a(4);
tbytes = a(5);
bbytes = a(6);
[a, count] = fread(id, 2, 'int16');
vers_id = a(1);
status = a(2);
[a, count] = fread(id, 1, 'int32');
nbheaders = a(1);

 BinaryStatus = fliplr(dec2bin(status));
 
 if ((BinaryStatus(3) == '0') & (BinaryStatus(4) == '0'))
    aux = 'int16'; 
 elseif ((BinaryStatus(3) == '1') & (BinaryStatus(4) == '0'))
    aux = 'int32';
 else
     aux = 'float';
 end
 
 if nargin == 2; nor =1; end
 
    for k=1:nblocks
        for m=1:nbheaders
          fread(id, nbheaders*14, 'int16');			%read in the block headers (nbheaders*28 bytes)
        end
      
        [b, count] = fread(id, ntraces*np, aux);  	%read in the actual data (ntraces*np)
        
          
        c1 = (b(1:2:end)); %detrend for varian is not good
        c2  =(b(2:2:end)); 
        struc{k}.fid = (c1 - i*c2)*exp(i*rp*pi/180);
        struc{k}.esp = fftshift(fft(struc{k}.fid,SI))*nor;
        struc{k}.sw = SW;
		struc{k}.tof = tof;
		struc{k}.lb = lb;
      end  
    
    
 
 
  
  