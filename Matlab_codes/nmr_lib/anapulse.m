function  U = anapulse(ang,tp,shape,freng,npf,ramp);
%function 
% analisa pulsos

shp = load([shape]);

np = length(shp);  dt = tp/np;   t = (1:np)*dt;

wr=pi*ang/(180* dt*sum(shp(:,2)) )*shp(:,2);

ph = shp(:,1);


if nargin == 6
dphidt = 2*pi*ramp;
ph = pi*ph/180 + dphidt*t';
end

w = 2*pi*(-freng:(2*freng/npf):freng);
 
for k=1:np
         Hrf{k} = + wr(k)*(cos(ph(k))*xgate/2 + sin(ph(k))*ygate/2);
end

roh = [1 0; 0 0];
  
for m=1:length(w)
     
     H = w(m)*zgate/2;

     U =1;
   
     for k=1:np
         
         U = expm(-i*(Hrf{k} + H)*dt)*U;
           
     end

      
      roh2 = U*roh*U';
      mz(m) = real(trace(zgate*roh2));
      mx(m) = trace((xgate)*roh2);
      my(m) = trace((ygate)*roh2);
      
      ff(m) = fidelmat(U,rotx(ang));
end


figure
plot((1e-3)*w/(2*pi),mz,'k',(1e-3)*w/(2*pi),mx,'r',(1e-3)*w/(2*pi),my,'b'); 
xlabel('Frequency (KHz)')
ylabel('Magnatization')
axis([-freng*1e-3 freng*1e-3 -1 1])

figure
plot((1e-3)*w/(2*pi),ff,'k'); 
axis([-freng*1e-3 freng*1e-3 -1 1])


