function alvo = func_rho_1q(rho_exp, chute, p)
%
%
% chi_2 = func_rho_2q(roh_exp, chute);
%
%

t = chute;

T_t = [t(1) 0; 
         t(3) + 1i * t(4)  t(2)];


rho_t = T_t' * T_t; rho_t = rho_t / trace(rho_t);


chi_2 = sum(sum(real(rho_t - rho_exp).^2 + imag(rho_t - rho_exp).^2)); 


if nargin > 2
    
    
    alvo = rho_t;
    
else
    
    alvo = chi_2;
    
end
    

