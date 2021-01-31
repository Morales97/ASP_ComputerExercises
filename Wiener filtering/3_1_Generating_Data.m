%% Generate white noise

K = 1000        % Num samples
sigma2_e = 1;
sigma2_w = 1;

e = sqrt(sigma2_e) * randn(K, 1);    % Kx1 column vector
w = sqrt(sigma2_w) * randn(K, 1);    % randn is gaussian distr. with zero mean


%% Generating a polynomial
% To define a polynomial A(q) we have 2 options:

% 1. Define directly its coefficients
% A(z) = z^2 - 4z + 3
A = [1 -4 3]    

% 2. Find the polynomial with the set of roots [p1 p2 ...]
% A(z) will be the same, 3 and 1 are the roots of A = [1 -4 3]  
A = poly([3 1])


%% Generate AR process
% The AR process x(n) is formed from white noise as
%   A(q)·x(n) = e(n)
% e.g., 
%   A(q) = 1 + a1·q^(-1)
%   which yields x(n) + a1·x(n-1) = e(n)
% To generate x(n) from white noise, find x(n) = 1/A(q)·e(n)

% 1. Define A(q). 
A = [1 0.5 0.8];    
% 2. Apply x(n) = 1/A(q)·e(n)
x = filter(1, A, e);

%% Generate colored noise

A_noise = [1 0.2];
v = filter(1, A_noise, e);

%% Generate noisy observations

y = x + v;

%% 
[SigmaYY, SigmaYx] = firw_cov_add(A, 1, A_noise, 1, 10);
