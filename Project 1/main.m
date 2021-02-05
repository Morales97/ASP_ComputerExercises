addpath(pwd)
cd ..
cd 'Spectral estimation'/
addpath(pwd)
cd ..
cd 'Wiener filtering'/
addpath(pwd)
cd ..
cd ..
cd mfiles/
addpath(pwd)
cd ..
[z,fs] = audioread('EQ2401project1data2021.wav');

%% FIR filter
[shatfir, thetahatfir] = p1_firw(z, z(1:8000), 40, 5);

%% Non-Causal filter
[shatnc, numnc, dennc] = p1_ncw(z, z(1:8000), 20, 10);

%% Causal filter
[shatc, numc, denc] = p1_cw(z, z(1:8000), 10, 10);

%% FIR AR filter
[shatfirAR, thetahatfirAR] = p1_firw_AR(z, z(1:8000), 10, 20, 15);

%% Non-Causal AR filter
[shatncAR, numncAR, denncAR] = p1_ncw_AR(z, z(1:8000), 30, 10);
