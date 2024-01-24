function struc = apodize(struc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(struc)
    t = ((0:length(struc{k}.fid)-1)')*(1/struc{k}.sw);
    SI = length(struc{k}.esp);
    struc{k}.fid = struc{k}.fid.*exp(-pi*struc{k}.lb.*t);
    struc{k}.esp = fftshift(fft(struc{k}.fid,SI));
end
   