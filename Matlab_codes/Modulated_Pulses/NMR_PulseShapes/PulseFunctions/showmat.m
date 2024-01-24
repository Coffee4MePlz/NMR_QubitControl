function showmat(fnun,roh);
%function showmat(fnum,roh);
%show matrix
%fnum = number of figure
%roh matrix to be plotted 

figure(fnun)

n = length(roh);
              
aux = [max(abs(real(roh))) max(abs(imag(roh)))];  
              
aux2 = [0 n+1 0 n+1 -max(aux) max(aux)];
              
subplot(121);

bar3(real(roh));

axis(aux2);

 texto1 = title('Real');
%texto1 = title('Parte Real','Color', 'w','fontsize',18);
 %'Case number # 3','Color', 'm'

axis off;
              
subplot(122);

bar3(imag(roh));

axis(aux2);

%texto1 = title('Parte Imaginária','Color', 'w','fontsize',18);
texto2 = title('Imaginary');

axis off;
