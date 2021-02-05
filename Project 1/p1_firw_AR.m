function [shat, thetahatfir] = p1_firw_AR(z, x, N, M_signal, M_noise)

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
    %[Ahat, sigma2hat] = ar_id(z(12000:30000), M_signal);
    %Ahat = Ahat';
    [Ahat, sigma2hat] = aryule(z, M_signal);
    
    % Get AR estimate for noise x(n)
    [Anoisehat, sigma2noisehat] = ar_id(x, M_noise);
    Anoisehat = Anoisehat';
    

    [SigmaZZ,SigmaZs] = firw_cov_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat, N);
    thetahatfir = SigmaZZ\SigmaZs;
    shat = filter(thetahatfir,1, z);
    
    % Plot estimated spectra
    
    %output
    [Ahatout, sigma2hatout] = aryule(shat, M_signal);

    w=linspace(0,pi);
    [mags,phases,ws]=dbode(1,Ahat,1,w);
    [magout,phases,wout]=dbode(1,Ahatout,1,w);
    [magv,phasev,wv]=dbode(1,Anoisehat,1,w);
    [magfir,phasefir,wfir]=dbode(thetahatfir',1,1,w);
    plt = semilogy(ws,mags.^2*sigma2hat,'b', wout,magout.^2*sigma2hatout,'r',wv,magv.^2*sigma2noisehat,':k', wfir, magfir.^2, '--');
    set(plt, 'LineWidth', 1.5)
    legend('Input','Output','Noise', 'FIR filter')
    title('Spectra')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
    
    