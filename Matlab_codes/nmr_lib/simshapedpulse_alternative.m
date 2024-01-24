function [U roh] = simshapedpulse(shape1,shape2,pw,phase1,phase2,ramp1,ramp2,power1,power2,roha)


phase1 = phase1*pi/180;

phase2 = phase2*pi/180;

global mol spectro;

     Ix1 = mol.Ich{1,1};   Iy1 = mol.Ich{1,2};
     Ix2 = mol.Ich{2,1};   Iy2 = mol.Ich{2,2};

shp = load([shape1]);

np = length(shp(:,3));  np2 = sum(shp(:,3));

dt = 1e-6*pw/np2;   t = ((1:np)*dt).*shp(:,3)';

wr1= 2*pi*power1*shp(:,2)/max(shp(:,2));
dphidt = 2*pi*ramp1;

ph1 = shp(:,1)*pi/180 +  phase1 + dphidt*t' ;

shp = load([shape2]);

np = length(shp(:,3));  np2 = sum(shp(:,3));

dt = 1e-6*pw/np2;   t = ((1:np)*dt).*shp(:,3)';

wr2= 2*pi*power2*shp(:,2)/max(shp(:,2));
dphidt = 2*pi*ramp2;

ph2 = shp(:,1)*pi/180 +  phase2 + dphidt*t' ;

U = 1;

for k=1:np 

         Hrf =   wr1(k)*(cos(ph1(k))*Ix1 + sin(ph1(k))*Iy1) + ...
                 wr2(k)*(cos(ph2(k))*Ix2 + sin(ph2(k))*Iy2);

         U = expm(-1i*(Hrf +  mol.Hint + mol.Hzee)*dt)*U;

end
[Ua roh]= delay(spectro.rof1*1e-6,roha); 
roh = U*roh*U'; 
[Ub roh]= delay(spectro.rof2*1e-6,roh);

U = Ub*U*Ua;
