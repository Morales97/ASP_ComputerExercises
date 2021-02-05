function [shat, numnc, dennc] = p1_ncw_AR(z, x, M_signal, M_noise)

% 
% [shat, thetahatfir] = p1_ncw_AR(z, x, M_signal, M_noise)
%
%   z           - Noisy sequence z(n) = s(n) + x(n)
%   x           - Noise samples x(n), when speaker is silent
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

    % Get spectra estimates
    [PhizsNum,PhizsDen,PhizzNum,PhizzDen] = ...
               spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
     
    % --- Compare Parametric vs Non-Parametric estimation ---  
    N = 20;
    % Blackman-Tuckey with Hamming window
    SigmaZzhat = xcovhat(z,z,N);
    PhizzNum_BT = [flip(SigmaZzhat(2:N)); SigmaZzhat];
    window = hamming(2*N-1);
    U = sum(window.^2)/length(window);
    PhizzNum_BT = (PhizzNum_BT .* window)'/U;
    PhizzDen_BT = [zeros(N-1,1); 1]';
    
    % Periodogram ensemble average (Bartlett's method)
    K = 135;
    interval = length(z) / K;
    SigmaZzhat_tmp = xcovhat(z(1:interval),z(1:interval),N);
    PhizzNum_BM = [flip(SigmaZzhat_tmp(2:N)); SigmaZzhat_tmp]';
    for i = 1:K-1
        t = z(interval*i:interval*(i+1));
        SigmaZzhat_tmp = xcovhat(t,t,N);
        PhizzNum_BM = PhizzNum_BM + [flip(SigmaZzhat_tmp(2:N)); SigmaZzhat_tmp]';
    end
    PhizzNum_BM = PhizzNum_BM / K;
    PhizzDen_BM = [zeros(N-1,1); 1]';
    
    w=linspace(0, pi, 512);
    [magp,~,wp]=dbode(PhizzNum,PhizzDen,1,w);
    [magbt,~,wbt]=dbode(PhizzNum_BT,PhizzDen_BT,1,w);
    [magbm,~,wbm]=dbode(PhizzNum_BM,PhizzDen_BM,1,w);
    figure
    plt = semilogy(wp, magp.^2, wbt, magbt.^2, wbm, magbm.^2);
    set(plt, 'LineWidth', 1.5)
    legend('Parametric estimation', 'Blackman-Tuckey (hamming)', 'Bartlett method')
    title('Z spectrum')
    % ---------------------------------------------------------
           
    % Use Blackman-Tuckey / Parametric estimation
    use_BT_spectrum = 1;
    use_BT_cross_spectrum = 0;
    if use_BT_spectrum
        PhizzNum = PhizzNum_BT;
        PhizzDen = PhizzDen_BT;
        SigmaXxhat = xcovhat(x,x,N);
        if use_BT_cross_spectrum
            PhixxNum = [flip(SigmaXxhat(2:N)); SigmaXxhat]';
            PhizsNum = PhizzNum - PhixxNum;
            PhizsDen = PhizzDen;
        end
    end
    
    % Get Non-causal filter
    numnc = conv(PhizsNum, PhizzDen);
    dennc = conv(PhizsDen, PhizzNum);

    shat = ncfilt(numnc, dennc, z);
    
    
    
    % --- PLOT ---
    [Aouthat, sigma2outhat] = aryule(shat, M_signal);
    
    w=linspace(0, pi, 512);
    [magz,~,wz]=dbode(1,Ahat,1,w);
    [mags,~,ws]=dbode(1,Aouthat,1,w);
    [magx,~,wx]=dbode(1,Anoisehat,1,w);
    [magnc,~,wnc]=dbode(numnc,dennc,1,w);
    figure
    plt = semilogy(wz,magz.^2*sigma2hat, 'b', ws, mags.^2*sigma2outhat, 'r', wx, magx.^2*sigma2noisehat, 'k:', wnc, magnc.^2, '--');
    set(plt, 'LineWidth', 1.5)
    title('Spectra')
    legend('Input z', 'Output shat', 'Noise', 'Non-causal freq response')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')

    
    
    