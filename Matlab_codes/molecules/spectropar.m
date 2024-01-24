%global - Declare variables as global
global spectro

spectro.chanel{1} = [1];
spectro.chanel{2} = [2];

% rof - delays (us)
spectro.rof1 = 10; %
spectro.rof2 = 24.875;
%SW - spectral width (Hz)
spectro.sw = 1996.8;
%Nfid - number of points of measurement
spectro.nfid  = 32*1024;
%Nesp - number of points to make FFT (fast fourier transform)
spectro.nesp = 32*1024;
%Spin observed - wich chanel is the observer
spectro.obs = 1;