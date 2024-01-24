function [rho_final,t_0] = find_rho_2q_v2(rho_exp, flag, target, figuras)
% Se figuras for diferente de zero o programa mostrara as figuras das
% matrizes antes e depois. 
% target --> estado esperado para chute inicial funciona junto com flag = 3
% flag - diferents cutes iniciais 
% flag = 1 chute inicial aleatório
% flag = 2 chute inicial ones normalizado
% flag = 3 chute inicial é o traget físico esperado

aaa = rho_exp;

fator = 10^7;

rho_med = (round(real(aaa) * fator) + 1i * round(imag(aaa) * fator)) / fator;

%rho_med; % mostra rho depois da primeira correcao


%------------------------------ finder ------------------------------------
%------ chute inical para construir o conjunto ----------------------------
%------ de parahmetros que entra no ajuste --------------------------------

if flag == 1
    t_0 = rand(1, 16);   % t's iniciais aleatorios
end

if flag == 2 
    t_0 = ones(1, 16);   % t's iniciais iguais a 1
end


if flag == 3  % ts iniciais igual ao target
    rr = target;
    % reescreve o target em termos uma uma matrix triagular inferior coef
    % t(j)
    t4 = sqrt(rr(4,4));
    if t4 == 0
        t9 = 0; t10 = 0; t13 =0; t14 = 0; t15 = 0; t16 = 0;
    else
        t9 = real(rr(4,3))/t4;
        t10 = imag(rr(4,3))/t4;
        t13 = real(rr(4,2))/t4;
        t14 = imag(rr(4,2))/t4;
        t15 = real(rr(4,1))/t4;
        t16 = imag(rr(4,1))/t4;
    end
    t3 = sqrt(rr(3,3) - (t9^2 + t10^2));
    if t3 == 0
        t7 = 0; t8 =0; t11 = 0; t12 = 0;
    else
        t7 = real((rr(3,2) - (t9 - t10*1i)*(t13 + t14*1i))/t3);
        t8 = imag((rr(3,2) - (t9 - t10*1i)*(t13 + t14*1i))/t3);
        t11 = real((rr(3,1) -(t9 - t10*1i)*(t15 + t16*1i))/t3);
        t12 = imag((rr(3,1) -(t9 - t10*1i)*(t15 + t16*1i))/t3);
    end
    t2 = sqrt(rr(2,2) - (t7^2 + t8^2 + t13^2 + t14^2));
    if t2 == 0
        t5 = 0; t6 = 0;
    else
        t5 = real((rr(2,1) - ((t7 - t8*1i)*(t11 + t12*1i) + (t13 - t14*1i)*(t15 + t16*1i)))/t2);
        t6 = imag((rr(2,1) - ((t7 - t8*1i)*(t11 + t12*1i) + (t13 - t14*1i)*(t15 + t16*1i)))/t2);
    end
    t1 = sqrt(rr(1,1) - (t5^2 + t6^2 + t11^2 + t12^2 + t15^2 + t16^2));
    t_0 = [t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16];
end

%------ mostra a matriz inicial que entra no ajuste -----------------------

rho_t = func_rho_2q(rho_med, t_0,1);   

%------ mostra o chi quadrado inicial -------------------------------------

chi_2_ini = func_rho_2q(rho_med, t_0);  


%------------- faz o ajuste - o resultado eh guardado em t_f --------------


options = optimset('MaxFunEvals', 2e11,'MaxIter', 1e11,'TolX', 1e-10);


t_f = fminsearch(@(t) func_rho_2q(rho_med, t), t_0, options);  % faz o ajuste


%------------- motra a matriz final ---------------------------------------


rho_f = func_rho_2q(rho_med, t_f, 1);  


%------------- mostra o chi quadrado final --------------------------------


chi_2_fin = func_rho_2q(rho_med, t_f);  % mostra o chi quadrado final


%----------------------------- final da funcao ----------------------------


rho_final = rho_f;


%------------- figura 1 --- matriz experimental e de entrada no ajuste ----


if figuras ~= 0


figure(998);

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


figure(999);

subplot(2,2,1);

bar3(real(rho_t)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;

subplot(2,2,2);

bar3(imag(rho_t)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;

subplot(2,2,3);

bar3(real(rho_f)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;

subplot(2,2,4);

bar3(imag(rho_f)); axis([0.5 4.5 0.5 4.5 -0.1 1.1]); axis off;


%--------------------------------------------------------------------------

end


