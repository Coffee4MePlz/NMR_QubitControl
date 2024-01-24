function sim_medida(rohtini,noise_size,T2,amplitudelim,nunfig,flag1);
% Simulacao da medida em RMN FID e Espectros com pulsos de peapracao
% 2 qubits
% Requer sim_fid.m e sim_fid2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sintax Simulacao dos Pulsos e Evolucao Livre
% rgpulse(pw1,ang1,phase1,roha)
% decrgpulse(pw2,ang2,phase2,roha)
% simpulse(pw1,pw2,ang1,ang2,phase1,phase2,roha)
% pw# tempo do pulso de 90, ang# angulo de rotacao, 
% phase# fasedo pulso (0 = X, 90 = Y, 180 = -X, 270 = -Y),
% rhoa estado antes do pulso
% delay((tempo/JJ),roha)
% tempo em unidades de 1/J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spectropar;
cloroformio;

% Definicoes uteis
id = [1 0; 0 1];
sz = [1 0; 0 -1]/2;
sx = [0 1; 1 0]/2;
sy = [0 -1i; 1i 0]/2;
ZZ = kron(sz,sz);
sminus = [0 1; 0 0];
splus = [0 0; 1 0];

% Atualizar com os dados do experiemento atual
pw1 = 11.8;   % Tempo pulso de 90 observer (microsegundo)
pw2 = 9.775;  % Tempo pulso de 90 decopler (microsegundo)
JJ = 215.10; % Freq. do acoplamento J em Hertz

te = 0;  % erro na rotacao
tff = 0; % erro na fase 

% Tomagrafia parte 1
for nr4 = 1:4
    % Pulso de preparacao da medida - Tomo parte 1
    roh = rohtini;
    if nr4 == 2; 
        [U roh] = rgpulse(pw1,90+te,90+tff,roh); % YI
    end
    if nr4 == 3
       [U roh] = rgpulse(pw1,90+te,0+tff,roh); % XI
    end
    if nr4 == 4
       [U roh] = simpulse(pw1,pw2,90+te,90+te,0+tff,0+tff,roh); % XX
    end   
    % Simula o Fid
    [fid, atempo] = sim_fid(roh,T2);
    afid1{nr4} = fid;
    atempo1{nr4} = atempo;
end

if flag1 == 1 % mostra ou não os fids
% Figuras para o Fid parte 1
figure(nunfig+1)
subplot(2,2,1)
plot(atempo1{1},real(afid1{1}))
t=title(['Part 1 - II']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);

subplot(2,2,2)
plot(atempo1{2},real(afid1{2}))
t=title(['Part 1 - YI']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);

subplot(2,2,3)
plot(atempo1{3},real(afid1{3}))
t=title(['Part 1 - XI']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);

subplot(2,2,4)
plot(atempo1{4},real(afid1{4}))
t=title(['Part 1 - XX']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);
end

% Tomagrafia parte 2 - medida no segundo qubit
for nr4 = 1:4
    % Pulso de preparacao da medida - Tomo parte 2
    roh = rohtini;
    if nr4 == 2; 
        [U roh] = rgpulse(pw1,90+te,90+tff,roh); % YI
    end
    if nr4 == 3
       [U roh] = rgpulse(pw1,90+te,0+tff,roh); % XI
    end
    if nr4 == 4
       [U roh] = simpulse(pw1,pw2,90+te,90+te,0+tff,0+tff,roh); % XX
    end     
    % Simula o Fid medindo no segundo qubit
    [fid, atempo] = sim_fid2(roh,T2);
    afid2{nr4} = fid;
    atempo2{nr4} = atempo;
end

if flag1 == 1 % mostra ou não os fids
% Figuras para o Fid parte 2
figure(nunfig+2)
subplot(2,2,1)
plot(atempo2{1},real(afid2{1}))
t=title(['Part 2 - II']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);

subplot(2,2,2)
plot(atempo2{2},real(afid2{2}))
t=title(['Part 2 - IY']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);

subplot(2,2,3)
plot(atempo2{3},real(afid2{3}))
t=title(['Part 2 - IX']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);

subplot(2,2,4)
plot(atempo2{4},real(afid2{4}))
t=title(['Part 2 - XX']);
t.FontSize = 12;
xlabel('t (s)', 'FontSize', 12);
ylabel('FID (a.u.)', 'FontSize', 12);
xlim([0 4]); ylim([-0.5 0.5]);
end

%% Simula os Espectros

% Espectros Parte I
for nr4 = 1:4;
    L = length(atempo1{nr4}); % numero de pontos de aquisicao do sinal do fid
    n = 2^nextpow2(L); % encontra a potencia de 2 mais proxima para FFT
    % Simulando ruido no fid
    afid1{nr4} = afid1{nr4} + noise_size*randn(size(fid)); % adiciona ruido
    aaux = afid1{nr4};
    Y = fft(aaux,n); % Fast Fourier Transform
    aY{nr4} = Y;
end

% Econtra os picos
% [pks,locs] = findpeaks(real(Y));

[tran1 tran1 ff1 opr] = transitions(mol,1);

af = ((ff1(1) - ff1(2))/7310)*(1:n) - ((ff1(1) - ff1(2))/7310)*n/2.0;
if flag1 == 0
    figure(nunfig+1)
else
    figure(nunfig+3) 
end
subplot(2,2,1)
plot(af,real(aY{1}));
t=title(['Part 1 - II']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)', 'FontSize', 12)
xlabel('Frequency (Hz)', 'FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])

subplot(2,2,2)
plot(af,real(aY{2}));
t=title(['Part 1 - YI']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)','FontSize', 12)
xlabel('Frequency (Hz)','FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])

subplot(2,2,3)
plot(af,real(aY{3}));
t=title(['Part 1 - XI']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)','FontSize', 12)
xlabel('Frequency (Hz)','FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])

subplot(2,2,4)
plot(af,real(aY{4}));
t=title(['Part 1 - XX']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)','FontSize', 12)
xlabel('Frequency (Hz)','FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])

% Espectros Parte II
for nr4 = 1:4;
    L = length(atempo2{nr4}); % numero de pontso de aquisicao do sinal do fid
    n = 2^nextpow2(L); % encontra a potencia de 2 mais proxima para FFT
    % Simulando ruido no fid o dobro de ruido para o carbono
    afid2{nr4} = afid2{nr4} + 2*noise_size*randn(size(fid)); % adiciona ruido
    aaux = afid2{nr4};
    Y2 = fft(aaux,n); % Fast Fourier Transform
    aY2{nr4} = Y2;
end

% Econtra os picos
% [pks,locs] = findpeaks(real(Y));

[tran1 tran1 ff1 opr] = transitions(mol,1);

af = ((ff1(1) - ff1(2))/7310)*(1:n) - ((ff1(1) - ff1(2))/7310)*n/2.0;

if flag1 == 0
    figure(nunfig+2)
else
    figure(nunfig+4) 
end

subplot(2,2,1)
plot(af,real(aY2{1}));
t=title(['Part 2 - II']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)', 'FontSize', 12)
xlabel('Frequency (Hz)', 'FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])

subplot(2,2,2)
plot(af,real(aY2{2}));
t=title(['Part 2 - YI']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)','FontSize', 12)
xlabel('Frequency (Hz)','FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])

subplot(2,2,3)
plot(af,real(aY2{3}));
t=title(['Part 2 - XI']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)','FontSize', 12)
xlabel('Frequency (Hz)','FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])

subplot(2,2,4)
plot(af,real(aY2{4}));
t=title(['Part 2 - XX']);
t.FontSize = 12;
ylabel('Amplitude (a.u.)','FontSize', 12)
xlabel('Frequency (Hz)','FontSize', 12)
set(gca,'XDir','reverse')
ylim([-amplitudelim amplitudelim])


