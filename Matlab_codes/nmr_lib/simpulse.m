function [U roh] = simpulse(pw1,pw2,ang1,ang2,phase1,phase2,roha)
% simulate simpulse (simutaneous square pulse on both obeserver and decoupler)
% [U roh] = simpulse(pw1,pw2,ang1,ang2,phase1,phase2,roha)
% pw1(2)  - > pulse width in microseconds for observer (decoupler)
% ang1(2) -> rotation angle for observer (decoupler)
% phase1(2) - > phase of the pulse in radiands for observer (decoupler)
% roha -> density matrix before the pulse
% roh -> density matrix after the pulse
% U ->  Pulse propagator (this incluies the free evoluiton induced by the rofs times)


global mol spectro;

phase1 = phase1*pi/180;
phase2 = phase2*pi/180;
ang1 = ang1*pi/180;
ang2 = ang2*pi/180;

         Ix1 = mol.Ich{1,1};   Iy1 = mol.Ich{1,2};
         Ix2 = mol.Ich{2,1};   Iy2 = mol.Ich{2,2};

if pw1 <= pw2; 

    Hrf1 =  ang1*(cos(phase1)*Ix1 + sin(phase1)*Iy1);
    Hrf2 =  ang2*(cos(phase2)*Ix2 + sin(phase2)*Iy2);
    tp = [(pw1-pw2)/2 pw2];
    aa = ((pw1-pw2)/2)/pw1; bb = pw2/pw1;
else 
    Hrf1 =  ang2*(cos(phase2)*Ix2 + sin(phase2)*Iy2);
    Hrf2 =  ang1*(cos(phase1)*Ix1 + sin(phase1)*Iy1);
    tp = [(pw2-pw1)/2 pw1];
    aa = ((pw2-pw1)/2)/pw2; bb = pw1/pw2;
end 
H = mol.Hzee + mol.Hint;
U = expm(-1i*(Hrf1*aa + (H)*tp(1)*1e-6)) * ...
     expm(-1i*(Hrf1*bb + Hrf2 + (H)*tp(2)*1e-6)) * ...
     expm(-1i*(Hrf1*aa +  (H)*tp(1)*1e-6));

[Ua roh]= delay(spectro.rof1*1e-6,roha); 
roh = U*roh*U'; 
[Ub roh]= delay(spectro.rof2*1e-6,roh);


 t = (tp(2) + 2*tp(1))*1e-6 + spectro.rof1*1e-6 + spectro.rof2*1e-6;
%  Td(1) = 0.29/4 ;
%  Td(2) = 0.16/4;
% 
% 
%  for mn = 1:1:2
%      
%  alfa(mn) = (1 + exp(-t/(2*Td(mn))))/2;
%  E1{mn} = [sqrt(alfa(mn)) 0; 0 sqrt(alfa(mn))];
%  E2{mn} = [sqrt(1 - alfa(mn)) 0; 0 -sqrt(1 - alfa(mn))];
% % 
% 
% end 
% Ef1 = kron(E1{1},E1{2});
% Ef2 = kron(E2{1},E2{2});
% roh = Ef1*roh*Ef1' + Ef2*roh*Ef2';

U = Ub*U*Ua;


