function [shat, thetahatfir] = p1_firw(z, x, N, M_signal, M_noise)

% 
% [shat, thetahatfir] = p1_firw_AR(z, x, N, ar_order)
%
%   z           - Noisy sequence z(n) = s(n) + x(n)
%   x           - Noise samples x(n), when speaker is silent
%   N           - Length of FIR filter
%   M_signal    - AR model order to estimate z(n)
%   M_noise     - AR model order to estimate noise
%
%   shat        - Estimate of s(n)
%   thetahatfir - Estimate of FIR coefficients
%
% Background noise cancellation with FIR filter of length N. The spectra of
% signal and noise are estimated with parametric methods, with an AR model
%
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get AR estimate for z(n)
    [Ahat, sigma2hat] = aryule(z, M_signal);
    
    % Get AR estimate for noise x(n)
    [Anoisehat, sigma2noisehat] = aryule(x, M_noise);
    
    % Autocorrelation of observations z(n)
    SigmaZz = ar2cov(Ahat, sigma2hat, N);
    SigmaZZ = toeplitz(SigmaZz);

    % Cross-correlation between z(n) and s(n)
    SigmaXx = ar2cov(Anoisehat, sigma2noisehat, N);
    SigmaSs = SigmaZz - SigmaXx;
    SigmaZs = SigmaSs;  
    
    % Compute FIR filter 
    thetahatfir = SigmaZZ\SigmaZs;
    
    % Estimate s(n)
    shat = filter(thetahatfir,1, z);
    
    
    % --- PLOTS ---
    [Ahatout, sigma2hatout] = aryule(shat, M_signal);

    w=linspace(0,pi);
    [mags,~,ws]=dbode(1,Ahat,1,w);
    [magout,~,wout]=dbode(1,Ahatout,1,w);
    [magv,~,wv]=dbode(1,Anoisehat,1,w);
    [magfir,~,wfir]=dbode(thetahatfir',1,1,w);
    plt = semilogy(ws, mags.^2*sigma2hat, 'b', ...
                   wout, magout.^2*sigma2hatout, 'r', ...
                   wv, magv.^2*sigma2noisehat, ':k', ...
                   wfir, magfir.^2, '--');
    set(plt, 'LineWidth', 1.5)
    legend('Input','Output','Noise', 'FIR filter')
    title('Spectra')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
    
    