%% AR model order selection
[z,fs] = audioread('EQ2401project1data2021.wav');
x = z(1:8000);      % Noise samples

M_lim = 400;
sig = z;
aic(sig, M_lim);

sig = x;
aic(sig, M_lim);