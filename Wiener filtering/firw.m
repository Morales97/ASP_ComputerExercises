
function [xhat, theta] = firw(y, SigmaYx, SigmaYY)
    % Generate optimal LMMSE estimator and return estimate - WORKS
    %
    % Input:
    %   y: sequence of observations of length K+1
    %   SigmaYx: cross-covariance between observations and process x 
    %   SigmaYY: covariance 

    theta = SigmaYY\SigmaYx;    % Equivalent to inv(SigmaYY) * SigmaYx 
    M = length(theta);
    y_extended = [zeros(M-1, 1); y];
    Gamma = toeplitz(y_extended); 
    Gamma = Gamma(M:end, 1:M);
    
    xhat = Gamma * theta;