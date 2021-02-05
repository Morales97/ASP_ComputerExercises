function [shat, numnc, dennc] = p1_ncw(z, x, N, M_noise)

% 
% [shat, thetahatfir] = p1_ncw_AR(z, x, M_signal, M_noise)
%
%   z           - Noisy sequence z(n) = s(n) + x(n)
%   x           - Noise samples x(n), when speaker is silent
%   N           - number of autocorrelation samples used to estimate spect.
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

    % Get noise spectrum
    [Anoisehat, sigma2noisehat] = aryule(x, M_noise);
    [PhixxNum, PhixxDen] = filtspec(1, Anoisehat, sigma2noisehat);
    
    w=linspace(0,pi, 512);
    [mag1,~,w1]=dbode(1,Anoisehat,1,w);
    [mag2,~,w2]=dbode(PhixxNum,PhixxDen,1,w);
    
    semilogy(w1, mag1.^2*sigma2noisehat, w2, mag2.^2)
    legend('AR', 'spec')
    figure
    % WARNING - for some reason, this plots are not the same (they should
    % be). Same curve, but differ by some factor  
    
    % Get Z spectrum
    SigmaZzhat = xcovhat(z,z,N);
    PhizzNum = [flip(SigmaZzhat(2:N)); SigmaZzhat]';
    PhizzDen = [zeros(N-1,1); 1]';
    [mag3,~,wv]=dbode(PhizzNum,PhizzDen,1,w);
    semilogy(wv, mag3)
    title('Z spectrum')
    figure


    % Estimate original signal spectrum
    SigmaXxhat = xcovhat(x,x,N);
    SigmaSshat = SigmaZzhat - SigmaXxhat;
    PhissNum = [flip(SigmaSshat(2:N)); SigmaSshat]';
    PhissDen = [zeros(N-1,1); 1]';
    [mag4,~,wv]=dbode(PhissNum,PhissDen,1,w);
    semilogy(wv, mag4)
    title('S spectrum')
    figure 
    
    numnc = conv(PhissNum, PhizzDen);
    dennc = conv(PhissDen, PhizzNum);

    
    shat = ncfilt(numnc, dennc, z);
    
    % plot
    [wbt, BT_spectrum_z] = BlackmanTuckey(z);
    [wbt, BT_spectrum_shat] = BlackmanTuckey(shat);
    
    w=linspace(0, pi, 512);
    [magv,phasev,wv]=dbode(1,Ahat',1,w);
    [magnc,~,wnc]=dbode(numnc,dennc,1,w);
    plt = semilogy(wbt, BT_spectrum_z(1:512).^2, wbt, BT_spectrum_shat(1:512).^2, wv, magv.^2*sigma2hat, wnc, magnc.^2, '--');
    set(plt, 'LineWidth', 1.5)
    title('Spectra')
    legend('Input z (BT)', 'Output shat (BT)','Noise', 'Non-causal freq response')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
