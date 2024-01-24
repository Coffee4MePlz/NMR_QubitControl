function f=fidelmat(r1,r2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Alexandre M. de Souza
% f=fidelmat(r1,r2) 
%
% Fidelity between r1 and r2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f =abs(trace(r1*r2'))/sqrt(trace(r1*r1')*trace(r2*r2'));
