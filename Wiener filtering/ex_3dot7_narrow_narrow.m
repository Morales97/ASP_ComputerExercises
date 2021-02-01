% Exercise 3.7 from Computer Exercises
% Narrow band signal + narrow band noise

%% White noise

K = 10000;        % Num samples
sigma2_e = 1;
sigma2_w = 1;

e = sqrt(sigma2_e) * randn(K, 1);    % Kx1 column vector
w = sqrt(sigma2_w) * randn(K, 1);    % randn is gaussian distr. with zero mean


%% Generate signal x
% Poles at mod=0.98 and nu=±1/4 (or e^(±pi/4))
module = 0.98;
angle = pi/4;
[pole_R, pole_I] = pol2cart(angle, module);     % Polar to cartesian
poles = [pole_R + 1i*pole_I, pole_R - 1i*pole_I];

% Generate polynomial A(q) with poles as its roots
A = poly(poles);
%zplane(roots(A))

% x = 1/A(q)·e(n)
x = filter(1, A, e);

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
M = 10;  % FIR filter order
N = M-1;
[SigmaYY, SigmaYx] = firw_cov_add(A, sigma2_e, Anoise, sigma2_w, N);
[PhixyNum, PhixyDen, PhiyyNum, PhiyyDen] = spec_add(A, sigma2_e, Anoise, sigma2_w);

[xhatnc, xhatc, xhatfir, numnc, dennc, numc, denc, thetahatfir] = est_add(x, v, N, A, sigma2_e, Anoise, sigma2_w, SigmaYx, SigmaYY);
figure()
spec_comp(A,sigma2_e,Anoise,sigma2_w,numnc,dennc,numc,denc,thetahatfir);

