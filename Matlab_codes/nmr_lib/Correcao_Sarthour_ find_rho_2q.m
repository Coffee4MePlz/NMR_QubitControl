function rho_final = find_rho_2q(rho_exp, figuras);
%
%
% rho_final = find_rho_2q(roh_exp, figuras);
%
%
% Se figuras for diferente de zero o programa mostrar? as figuras das
% matrizes antes e depois. 
%
%

f = min(diag(rho_exp));


if f > 0
    
    f = 0;
    
end;


aaa = rho_exp - f * eye(4); 


aaa = aaa / trace(aaa);


fator = 10^7;


rho_med = (round(real(aaa) * fator) + 1i * round(imag(aaa) * fator)) / fator;


rho_med % mostra rho depois da primeira correcao


%------------------------------ finder ------------------------------------


%------ chute inical para construir o conjunto ----------------------------

%------ de parahmetros que entra no ajuste --------------------------------


t_0 = rand(1, 16);   % t's iniciais aleatorios


t_0 = ones(1, 16);   % t's iniciais iguais a 1


%------ mostra a matriz inicial que entra no ajuste -----------------------


rho_t = func_rho_2q(rho_med, t_0, 1)   


%------ mostra o chi quadrado inicial -------------------------------------


chi_2_ini = func_rho_2q(rho_med, t_0)  


%------------- faz o ajuste - o resultado eh guardado em t_f --------------


options = optimset('MaxFunEvals', 2e11,'MaxIter', 1e11,'TolX', 1e-10);


t_f = fminsearch(@(t) func_rho_2q(rho_med, t), t_0, options);  % faz o ajuste


%------------- motra a matriz final ---------------------------------------


rho_f = func_rho_2q(rho_med, t_f, 1)  


%------------- mostra o chi quadrado final --------------------------------


chi_2_fin = func_rho_2q(rho_med, t_f)  % mostra o chi quadrado final


%----------------------------- final da funcao ----------------------------


rho_final = rho_f;


%------------- figura 1 --- matriz experimental e de entrada no ajuste ----


if figuras ~= 0


figure(1);

subplot(2,2,1);


escala = max(diag(rho_exp)); escala = escala + escala * 0.1;


bar3(real(rho_exp)); axis([0.5 4.5 0.5 4.5 -escala escala]); axis off;

subplot(2,2,2);

bar3(imag(rho_exp)); axis([0.5 4.5 0.5 4.5 -escala escala]); axis off;

subplot(2,2,3);

bar3(real(rho_med)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;

subplot(2,2,4);

bar3(imag(rho_med)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;



%------------- figura 2 --- matriz inicial (chute) e final do ajuste ------


figure(2);

subplot(2,2,1);

bar3(real(rho_t)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;

subplot(2,2,2);

bar3(imag(rho_t)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;

subplot(2,2,3);

bar3(real(rho_f)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;

subplot(2,2,4);

bar3(imag(rho_f)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;


%--------------------------------------------------------------------------

end;


