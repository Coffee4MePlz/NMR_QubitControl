function plotfit(fignum,spec1,frr)
% show fits of spectra 
% fignum - > number of figure
% spec1 - > structure with the sectrum 
% frr = [fr1 fr2] - > the figure will show tha range of of 
%freqeuncy between fr1 and fr2 


dt = 1/spec1.sw;

SI = length(spec1.esp);

fr1 = (1/dt)*(0:SI-1)/SI - spec1.sw/2 ;
fr2 = min(fr1):0.01:max(fr1);

freq = spec1.fr;
for k=1:length(freq)
    funlist{2*k-1} = @(c,x) (1/(2*pi*c(1)))./((x - freq(k)).^2 + (1/(2*pi*c(1)))^2);
    funlist{2*k}   = @(c,x) (x-freq(k))./((x - freq(k)).^2 + (1/(2*pi*c(1)))^2);
end

ylor = 0;
for k=1:length(freq)
    
    ylor =  ylor + spec1.ILP(2*k-1)*funlist{2*k-1}(spec1.INLP(1),fr2) + ...
            spec1.ILP(2*k)*funlist{2*k}(spec1.INLP(1),fr2);
end 
ylor = ylor+spec1.ILP(2*k+1);
figure(fignum)
plot(fr1,real(spec1.esp),'k',fr2,real(ylor),'r')
if nargin == 3 ; xlim(frr); end
ylabel('AMLITUDE (A.U.)')
xlabel('FREQUENCY (Hz)')
set(gca,'XDir','reverse')


