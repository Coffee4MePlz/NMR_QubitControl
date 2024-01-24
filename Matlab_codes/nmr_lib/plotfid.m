function plotfid(fignum,spec1)
% show FID
% fignum - > number of figure
% spec1 - > structure with the fid

t1 = (0:length(spec1.fid)-1)*(1/spec1.sw);

figure(fignum)

plot(t1,real(spec1.fid),'b')
hold on
plot(t1,imag(spec1.fid),'r')

ylabel('AMLITUDE (A.U.)')
xlabel('TEMPO (s)')
