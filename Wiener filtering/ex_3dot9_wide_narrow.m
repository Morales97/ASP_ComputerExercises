% Exercise 3.9 from Computer Exercises
% WIDE band signal + NARROW band noise

%% White noise

K = 10000;        % Num samples
sigma2_e = 1;
sigma2_w = 1;

e = sqrt(sigma2_e) * randn(K, 1);    % Kx1 column vector
w = sqrt(sigma2_w) * randn(K, 1);    % randn is gaussian distr. with zero mean

%% Let signal x be white
A = 1;
x = e;

%% Generate noise v
% Poles at mod=0.98 and nu=±1/3
module = 0.98;
angle = pi/3;
[pole_R, pole_I] = pol2cart(angle, module);     
poles = [pole_R + 1i*pole_I, pole_R - 1i*pole_I];

% Generate polynomial Anoise(q) with poles as its roots
Anoise = poly(poles);

% v = 1/Anoise(q)·w(n)
v = filter(1, Anoise, w);

%% Generate noisy observations
y = x + v;

%% Compute correlation and spectra, filters, estimates, and compare them
M = 3;  % FIR filter order
N = M-1;
[SigmaYY, SigmaYx] = firw_cov_add(A, sigma2_e, Anoise, sigma2_w, N);
[PhixyNum, PhixyDen, PhiyyNum, PhiyyDen] = spec_add(A, sigma2_e, Anoise, sigma2_w);

[xhatnc, xhatc, xhatfir, numnc, dennc, numc, denc, thetahatfir] = est_add(x, v, N, A, sigma2_e, Anoise, sigma2_w, SigmaYx, SigmaYY);
figure()
spec_comp(A,sigma2_e,Anoise,sigma2_w,numnc,dennc,numc,denc,thetahatfir);

%% Plot impulse responses
x_max = 20;
t = [0:x_max]';

% Causal filter impulse response
% Option 1: built-in function
% impz(numc,denc,x_max)

% Option 2: generating a delta and filtering it
figure
s = (t==0);
causal_resp = filter(numc, denc, s);
stem(t, causal_resp, 'filled')
hold on 

% FIR filter impulse response
s = [thetahatfir; zeros(x_max - length(thetahatfir) + 1, 1)];
stem(t, s, 'filled')
leg_fir = sprintf('FIR (%d)', M);
legend('Causal', leg_fir)
legend
title('Filter Impulse response')

% Non-Causal impulse response
figure()
t = [-x_max:x_max]';
s = (t==0);
nc_resp = ncfilt(numnc, dennc, s);
stem(t, nc_resp, 'filled')
title('Non-causal impulse response')

