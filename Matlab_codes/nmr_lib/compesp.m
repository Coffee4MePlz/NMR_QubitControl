function compesp(fignum,spec1,spec2,frr)
% compare two spectra
% fignum - > number of figure
% spec1 - > structure with the sectrum 1
% spec2 - > structure with the sectrum 1
% frr = [fr1 fr2] - > the figure will show tha range of of 
%freqeuncy between fr1 and fr2 

dt = 1/spec1.sw;

SI = length(spec1.esp);

fr1 = (1/dt)*(0:SI-1)/SI - spec1.sw/2 ;

dt = 1/spec2.sw;

SI = length(spec2.esp);

fr2 = (1/dt)*(0:SI-1)/SI - spec2.sw/2 ;

figure(fignum)
plot(fr1,real(spec1.esp),'b');
hold on
plot(fr2,real(spec2.esp),'r');
hold off
if nargin == 4 ; xlim(frr); end
ylabel('AMLITUDE (A.U.)')
xlabel('FREQUENCY (Hz)')
set(gca,'XDir','reverse')


