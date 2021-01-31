
function [SigmaYY, SigmaYx] = firw_cov_add(A, sigma2, Anoise, sigma2noise, N)
    % Find covariance matrix of Y and cross-covariance vector Yx - WORKS
    %
    % Model:
    %   x is an AR process defined by A
    %   v is AR noise defined by Anoise
    %   y = x + v
    %
    % Input:
    %   A: AR polynomial
    %   sigma2: variance of white noise used to generate x
    %   Anoise: AR polynomial of noise
    %   sigma2noise: variance of white noise used to generate v
    %   N: length of Y = [y(n) ... y(n-N+1)]
    %
    % Output:
    %   SigmaYY: covariance matrix of Y                 (N+1, N+1)
    %   SigmaYx: cross-covariance between Y and x(n)    (N+1, 1)


    % Since x and v are uncorrelated, covariance of y is cov of x + cov of v
    cov_x = ar2cov(A, sigma2, N);
    cov_v = ar2cov(Anoise, sigma2noise, N);
    cov_y = cov_x + cov_v;
    
    % Covariance matrix of Y has a toeplitz form
    SigmaYY = toeplitz(cov_y);
    
    % Cross-covariance vector Yx is the covariance of x (x, v independent)
    SigmaYx = cov_x;
    
    