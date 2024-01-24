function [U roh] = delay_Par(tp,roha,L_Mol)
%simualte a delay
%[U roh] = delay(tp,roha)
% tp time in seconds
% roha -> density matrix before the delay
% roh -> density matrix after the delay
% U ->  Propagator of the delay

U = expm(-1i*(L_Mol(:,:,2) + L_Mol(:,:,1))*tp);

roh = U*roha*U';

% % end
% 
% 
%  t = tp;
%  Td(1) = 1008.9 ;
%  Td(2) = 5000.9;
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
 end