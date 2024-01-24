function plotesp(fignum,spec1,frr)
% show spectrum
% fignum - > number of figure
% spec1 - > structure with the sectrum
% frr = [fr1 fr2] - > the figure will show tha range of of 
%freqeuncy between fr1 and fr2 

dt = 1/spec1.sw;

SI = length(spec1.esp);

fr1 = (1/dt)*(0:SI-1)/SI - spec1.sw/2 ;


figure(fignum)
plot(fr1,real(spec1.esp),'k')
if nargin == 3 ; xlim(frr); end
ylabel('AMLITUDE (A.U.)')
xlabel('FREQUENCY (Hz)')
set(gca,'XDir','reverse')

