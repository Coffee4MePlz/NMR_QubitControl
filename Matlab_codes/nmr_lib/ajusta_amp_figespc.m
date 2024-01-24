function ajusta_amp_figespc(amplitudelim,nunfig,flag1);
%% Ajuste de amplitude do espectro - executar essa seção para ajustar a visualizacao (opicional)

% amplitudelim  ------ limite da amplitude para o grafico do espectro (u.a.)
%                     % analogo ao vscale do VNMRJ - ajustar para visualizar
%                     % melhor
if flag1 == 0
    figure(nunfig+1)
else
    figure(nunfig+3) 
end
subplot(2,2,1); ylim([-amplitudelim amplitudelim])

subplot(2,2,2); ylim([-amplitudelim amplitudelim]);

subplot(2,2,3); ylim([-amplitudelim amplitudelim]);

subplot(2,2,4); ylim([-amplitudelim amplitudelim]);

% Espectros Parte II
if flag1 == 0
    figure(nunfig+2)
else
    figure(nunfig+4) 
end
subplot(2,2,1); ylim([-amplitudelim amplitudelim]);

subplot(2,2,2); ylim([-amplitudelim amplitudelim]);

subplot(2,2,3); ylim([-amplitudelim amplitudelim]);

subplot(2,2,4); ylim([-amplitudelim amplitudelim]);
% 