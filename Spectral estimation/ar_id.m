function [Ahat, sigma2hat] = ar_id(y, N)

% [Ahat,sigma2hat]=ar_id(y,N)
%
%	y			- Data sequence
%	N			- Model order 
%	Ahat			- Estimate of AR polynomial [1, a1, ..., aN]
%	sigma2hat		- Estimate of noise variance 
%
%
%  ar_id: Identification of AR model
%
%         Model: y(n)+a_{1}y(n-1)+...+a_{N}y(n-N)=e(n)
%
%     
%     Author: Daniel Morales
%
% NOTE: This function is equivalent to the built-in function for Yule-Walker
% [Ahat, sigma2hat] = aryule(z,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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