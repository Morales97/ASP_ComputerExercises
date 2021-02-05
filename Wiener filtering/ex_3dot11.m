
clear
load signal1.mat

M = 20;  % FIR filter order
N = M-1;
[SigmaYY, SigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noise, N);
[xhatnc, xhatc, xhatfir, numnc, dennc, numc, denc, thetahatfir] = est_add(x, v, N, A, sigma2, Anoise, sigma2noise, SigmaYx, SigmaYY);

%sound(y, fs)
%sound(xhatc,fs)
%sound(xhatnc,fs)
%sound(xhatfir,fs)