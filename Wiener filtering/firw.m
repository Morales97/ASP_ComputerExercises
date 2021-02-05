function [xhat,theta] = firw(y, SigmaYx, SigmaYY)

%
% [xhat,theta] = firw(y, SigmaYx, SigmaYY)
%	
%	y       - y(n)=x(n)+v(n)
% 	SigmaYY	- E[Y(n) (Y(n))']
%	SigmaYx	- E[Y(n) x(n)]
%	
% 	xhat	- FIR Wiener estimate of x(n) from y(n)
% 	theta	- FIR Wiener filter.
%	
%
%  firw: FIR Wiener estimate of x(n) from y(n)
%     
%     
%     Author: Daniel Morales
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    theta = SigmaYY\SigmaYx;    % Equivalent to inv(SigmaYY) * SigmaYx 
    M = length(theta);
    y_extended = [zeros(M-1, 1); y];
    Gamma = toeplitz(y_extended); 
    Gamma = Gamma(M:end, 1:M);
    
    xhat = Gamma * theta;