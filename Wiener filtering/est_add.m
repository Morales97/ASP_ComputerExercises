function [xhatnc,xhatc,xhatfir,numnc,dennc,numc,denc,thetahatfir] =...
    est_add(x, v, N, Ahat, sigma2hat, Anoisehat,...
    sigma2noisehat,SigmaYxhat, SigmaYYhat)
%
% [xhatnc,xhatc,xhatfir,numnc,dennc,numc,denc,thetahatfir] = 
%     est_add(x, v, N, Ahat, sigma2hat, Anoisehat,
%       sigma2noisehat,SigmaYxhat, SigmaYYhat)
%	
%	x			 - AR Signal 
%	v			 - AR Noise, y(n)=x(n)+v(n)
%	N			 - Length of the FIR Wiener filter
%	Ahat,sigma2hat 		 - Estimated or true parameters of x
%	Anoisehat,sigma2noisehat - Estimated or true parameters of v
%	SigmaYxhat		 - E[Y(n) x(n)]
% 	SigmaYYhat		 - E[Y(n) (Y(n))']
%	
% 	xhatnc		- Non-causal Wiener estimate of x
% 	xhatc		- Causal Wiener estimate of x
% 	xhatfir		- FIR Wiener estimate of x
% 	numnc,dennc	- Non-causal Wiener filter
% 	numc,denc	- Causal Wiener filter
% 	thetahatfir	- FIR Wiener filter
%	
%
%  est_add: Estimate using the three different Wiener filters.
%     Plot the 30 first samples of each estimate, x and y.
%     Calculate the MSE of the estimates normalized with the
%     inverse of the noise variance.
%     
%     Author: Daniel Morales
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Noisy observations
    y = x + v;
     
    % Get power spectrum of YY and XY
    [PhixyNum, PhixyDen, PhiyyNum, PhiyyDen] = spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
    
    % FIR Wiener
    [xhatfir, thetahatfir] = firw(y, SigmaYxhat, SigmaYYhat);
    
    % Non-causal Wiener
    [xhatnc, numnc, dennc] = ncw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);
    
    % Causal Wiener
    [xhatc, numc, denc] = cw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen, 0);     % Not sure if the prediction horizon m should correspond to parameter N, or if it should be set to 0. Then, what's N for?
    
    % Calculate MSE for the estimates, normalized by variance of noise v
    MSE_fir = mean((xhatfir - x).^2) / sigma2noisehat;
    MSE_nc = mean((xhatnc - x).^2) / sigma2noisehat;
    MSE_c = mean((xhatc - x).^2) / sigma2noisehat;
    
    leg_fir = sprintf('FIR est. (MSE = %.3f)', MSE_fir);
    leg_nc = sprintf('Non-causal est. (MSE = %.3f)', MSE_nc);
    leg_c = sprintf('Causal est. (MSE = %.3f)', MSE_c);

    
    figure()
    hold on
    x_axis = 1:1:30;
    
    plot(x_axis, x(1:30), 'DisplayName', 'Signal x')
    plot(x_axis, xhatfir(1:30), '-o', 'DisplayName', leg_fir)
    plot(x_axis, xhatnc(1:30), '-.', 'DisplayName', leg_nc)
    plot(x_axis, xhatc(1:30), '--', 'DisplayName', leg_c)
    plot(x_axis, y(1:30), ':', 'DisplayName', 'Observations y')

    legend