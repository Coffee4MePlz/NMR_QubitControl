function alvo = func_rho_2q(rho_exp, chute, p);
%
%
% chi_2 = func_rho_2q(roh_exp, chute);
%
%

t = chute;

T_t = [t(1) 0 0 0; 
         t(5) + 1i * t(6) t(2) 0 0; 
         t(11) + 1i * t(12) t(7) + 1i * t(8)  t(3) 0; 
         t(15) + 1i * t(16) t(13) + 1i * t(14) t(9) + 1i * t(10) t(4)];


rho_t = T_t' * T_t; rho_t = rho_t / trace(rho_t);


chi_2 = sum(sum(real(rho_t - rho_exp).^2 + imag(rho_t - rho_exp).^2)); 


if nargin > 2
    
    
    alvo = rho_t;
    
else
    
    alvo = chi_2;
    
end;
    

