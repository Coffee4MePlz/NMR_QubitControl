% Open the file
fileID = fopen('/home/user/nmr/Matlab_codes/Modulated_Pulses/NMR_PulseShapes/NMR_PulseShapes_output/last1simFSfullModRY90IDSA1Obs.RF', 'r');
%fileID = fopen('/home/user/nmr/Matlab_codes/Modulated_Pulses/NMR_PulseShapes/NMR_PulseShapes_output/gaussiansRFs/singlePhgaussModRY90IDpw200np200Obs.RF', 'r');
% Skip commented lines
line = fgets(fileID);
while line(1) == '#'
    line = fgets(fileID);
end

% Read the data
data = textscan(fileID, '%f %f %f');

% Close the file
fclose(fileID);

% Extract phase, amplitude, and timesteps
phase = data{1};
amplitude = data{2};
timesteps = data{3}; % You may or may not use this, depending on your needs
t = linspace(0,sum(timesteps)*1e-6, length(timesteps)+1);
amplitude = (heaviside(t - 2*1e-5) - heaviside(t - 3.1*1e-5));
phase = ones(1,length(t))*90;
phase(1) = 360;


%amplitude_2 = (heaviside(t - 2*1e-5) - heaviside(t - 5.5*1e-5))*5;

plotPulseFFT(amplitude, [], t,100,1,molfreq) %, numfig2) 
%plotPulseFFT(amplitude_2, [], t,10,2,molfreq) %, numfig2) 
disp('done')



