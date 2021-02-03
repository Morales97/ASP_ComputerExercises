function [Ahat, sigma2hat] = ar_id(y, N)
    % Find parameters for the AR-N model that generates y using one-step
    % ahead prediction
    % This function is equivalent to the built-in function for Yule-Walker
    % [Ahat, sigma2hat] = aryule(z,N)
    %
    % Model
    %   we assume y to be generated as an AR process of order N
    %   y(n) + a1·y(n-1) + ... + aN·y(n-N) = e(n)
    %   then, we can obtain AR parameters as the coef. of a FIR Wiener if
    %   we define x(n) = y(n+1)
    %   xhat(n) = A(q)·Y(n)
    %   yhat(n+1) = a1·y(n) + a2·y(n-1) + ... + aN·y(n-N+1)
    %
    % Input
    %   y: sequence, suposedly an AR-N process
    %   N: order of the AR model
    %
    % Output
    %   Ahat: AR parameters estimates, of the form [1 a1 a2 ... aN] 
    %   sigma2hat: driving white noise variance estimate
    
    % find covariance of Y
    SigmaYYhat = covhat(y, N);
    
    % find cross-covariance between Y and x
    % this is equivalent to covariance of Y shifted by 1 sample
    M = length(y);
    x = y(2:M);
    SigmaYxhat = xcovhat(x,y(1:M-1),N);
    
    % Find FIR coefficients
    [yhat, theta] = firw(y, SigmaYxhat, SigmaYYhat);
    Ahat = [1; -theta];
    
    % Estimate variance of driving noise
    sigma2hat = sum((y(2:M)-yhat(1:M-1)).^2) / (M-1);