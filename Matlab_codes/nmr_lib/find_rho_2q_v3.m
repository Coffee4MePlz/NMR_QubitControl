function [rho_final,fmaior,t_0_melhor] = find_rho_2q_v3(rho_exp, n_trials, target, figuras)
% Se figuras for diferente de zero o programa mostrara as figuras das
% matrizes antes e depois. 
% Repete n_trials vezes o find_rho_2q_v2 com chute inicial aleatorio e
% retorna a melhor fidelidade com realacao ao estado target
%

fmaior = 0.0;
for jj = 1:n_trials
    [rho_ml,t_0] = find_rho_2q_v2(rho_exp, 1, target, figuras);
    fiddml =  fidelmat(rho_ml,target);
    %fiddml = Fidelity(target,rho_ml); % Uhlmann fidelity
    if jj == 1
        %fmaior = fiddml;
        rho_ml_melhor = rho_ml;
        %t_0_melhor = t_0;    
    end
    if fiddml > fmaior
        fmaior = fiddml
        rho_ml_melhor = rho_ml;
        t_0_melhor = t_0;
    end
end

rho_final = rho_ml_melhor;

end