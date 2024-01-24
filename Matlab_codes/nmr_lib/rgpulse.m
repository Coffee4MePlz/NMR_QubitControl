function  [U roh] = rgpulse(pw,ang,phase,roha);
% simulate rgpulse (square pulse on obeserver)
%[U roh] = rgpulse(pw,ang,phase,roha);
% pw - > pulse width in microseconds
% ang -> rotation angle degree
% phase - > phase of the pulse in degree
% roha -> density matrix before the pulse
% roh -> density matrix after the pulse
% U ->  Pulse propagator (this incluies the free evoluiton induced by the rofs times)

global mol spectro;

ang = ang*pi/180;
phase = phase*pi/180;

Hrf =  (cos(phase)*mol.Ich{1,1} + sin(phase)*mol.Ich{1,2});

U = expm(-1i*(ang*Hrf + (mol.Hint + mol.Hzee)*pw*1e-6));

[Ua roh]= delay(spectro.rof1*1e-6,roha); 
roh = U*roh*U'; 
[Ub roh]= delay(spectro.rof2*1e-6,roh);


%  t = pw *1e-6 + spectro.rof1*1e-6 + spectro.rof2*1e-6;
%  Td(1) = 3.2 ;
%  Td(2) = 3.2;
% 
% 
% % for mn = 1:1:3
% %     
% % alfa(mn) = (1 + exp(-t/(2*Td(mn))))/2;
% % E1{mn} = [sqrt(alfa(mn)) 0; 0 sqrt(alfa(mn))];
% % E2{mn} = [sqrt(1 - alfa(mn)) 0; 0 -sqrt(1 - alfa(mn))];
% % 
% % end 
% % Ef1 = kron(E1{1},kron(E1{2},E1{3}));
% % Ef2 = kron(E2{1},kron(E2{2},E2{3}));
% % roh = Ef1*roh*Ef1' + Ef2*roh*Ef2';
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
 

