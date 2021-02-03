
[z,fs] = audioread('EQ2401project1data2021.wav');

%% FIR filter
[shatfir, thetahatfir] = p1_firw(z, z(1:8000), 50, 10);

%% Non-Causal filter
[shatnc, numnc, dennc] = p1_ncw(z, z(1:8000), 20, 10);