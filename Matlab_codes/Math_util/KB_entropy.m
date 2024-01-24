
%%
clc;
clear all;
close all

rho = [1.0 0.0;0.0 0.0]
M{1} = [1.0 0.0;0.0 0.0];
M{2} = [0.0 1.0; 0.0 0.0];
rho1 = M{1}*rho*M{1}'+M{2}*rho*M{2}'
sigma = [0.5 0.0 ; 0.0 0.5];

rel_entropy(rho,sigma)

%%_____Fucntions_______%%

%%_____Entropia von Neuman_______%%
function ent = vn_entropy(rho);
lambda = round(eig(rho),6);
for i = 1:length(lambda)
    if lambda(i) < 0
        disp('Error!!! neg. eigenvals - vN entropy')
        return
    end
end
% S(rho) = -Tr[rho ln(rho)] = -Sum_j (lambda_j ln lambda_j)
S = 0.0;
for i = 1:length(lambda)
    if lambda(i) > 0
        S = S - lambda(i)*log(lambda(i));
    end    
end   
ent = S;
end

%%_____Entropia relativa (Kullback)_______%%
function ent = rel_entropy(rho,sigma);
% S(rho/sigma) = Tr[rho*ln(rho)]-Tr[rho*ln(sigma)]
% 1º Dado rho e sigma (operadores), calcule os respectivos autovalores
% autovalores de rho
lambda_rho = round(eig(rho),6); 
for i = 1:length(lambda_rho)
    if lambda_rho(i) < 0
        disp('Error!!! neg. eigenvals - vN entropy')
        return
    end
end
%autovalores de sigma
lambda_sigma = round(eig(sigma),6); 
for i = 1:length(lambda_sigma)
    if lambda_sigma(i) < 0
        disp('Error!!! neg. eigenvals - vN entropy')
        return
    end
end
% entropia de von Neuman de rho
% S(rho) = -Tr[rho ln(rho)] = -Sum_j (lambda_j ln lambda_j)
Svn = 0.0;
for i = 1:length(lambda_rho)
    if lambda_rho(i) > 0
        Svn = Svn - lambda_rho(i)*log(lambda_rho(i));
    end    
end   

% relative trace = Tr[rho*ln(sigma)] = Sum_i (lam_rho*ln(lam_sigma))
rel_trace = 0.0;
for k =1:length(lambda_sigma)
    if lambda_sigma(k) > 0
        rel_trace = rel_trace + lambda_rho(k)*log(lambda_sigma(k));
    end
% entropia relativa S(rho/sigma)=-Sv(rho)-rel_trace
    
end

ent = -Svn - rel_trace;


end



