function [fid, atempo] = sim_fid(roh,T2);
% function [fid, atempo] = sim_fid(roh,T2)
% Simula o fid com decaimento exponecial 1/T2
% fid -> intencidade do sinal medido
% atempo -> pontos de aquisicao no dominio do tempo
% sempre medindo no qubit A no tx

%Simula o fid
sminus = [0 1; 0 0]; % sigma_{-} |0><1| observavel a ser medido
rohtini = roh; % inicializa o estado para simular fid
tempo = 0.0; % tempo inicial coleta do fid
sw = 1996; % taxa de aquisicao do Fid - spectral window in Hz (espaço das freq.)
dt = 1.0/sw; % incremento no tempo 
Np = 2^13; % numero de pontos para aquisicao 2^13 = 8195 ou 2^12 = 4096
G=1/T2;   % taxa de decaimento
for i = 1:Np
    roh = rohtini; %inicializa o estado
    [U roh] = delay(tempo,roh); % simula a evolucao livre no t
    aux = PartialTrace(roh,2); % estado reduzido qubit A
    fid(i) = trace(aux*sminus')*exp(-G*tempo); % fid observavel medido no qubit(a) com atenuação do sinal
    atempo(i) = tempo;  % quarda  o tempo
    tempo = tempo + dt; % incremento no tempo
end
%end of function [fid, atempo] = sim_fid(roh,T2)
