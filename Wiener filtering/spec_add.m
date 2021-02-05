function [PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
               spec_add(A, sigma2, Anoise, sigma2noise)

%
% [PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
%                 spec_add(A, sigma2, Anoise, sigma2noise)
%	
%	A           - AR model for the signal x(n), A(q)x(n)=w(n)
%	sigma2		- E[w(n)*w(n)]
%	Anoise		- AR model for the noise v(n), Anoise(q)v(n)=e(n)
%	sigma2noise	- E[e(n)*e(n)]
%	
% 	PhixyNum,PhixyDen	- Cross-spectrum between x(n) and y(n)
% 	PhiyyNum,PhiyyDen	- Spectrum of y(n)
%	
%  spec_add: Calculate spectrum and cross-spectrum for y(n)=x(n)+v(n)
%     
%     
%     Author: Daniel Morales
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    [PhixxNum, PhixxDen] = filtspec(1, A, sigma2);
    [PhivvNum, PhivvDen] = filtspec(1, Anoise, sigma2noise);

    % x and v uncorrelated, y is addition of their sectra
    [PhiyyNum, PhiyyDen] = add(PhixxNum, PhixxDen, PhivvNum, PhivvDen);
    
    % x and v uncorrelated, cross-spectrum xy is spectrum of x
    PhixyNum = PhixxNum;
    PhixyDen = PhixxDen;