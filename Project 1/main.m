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
cd data/
addpath(pwd)
cd ..
[z,fs] = audioread('EQ2401project1data2021.wav');

%% FIR filter
[shatfir, thetahatfir] = p1_firw(z, z(1:8000), 20, 10);

%% Non-Causal filter
[shatnc, numnc, dennc] = p1_ncw(z, z(1:8000), 20, 10);

%% Causal filter
[shatc, numc, denc] = p1_cw(z, z(1:8000), 10, 10);