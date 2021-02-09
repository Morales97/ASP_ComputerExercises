function [shat, numnc, dennc] = p1_ncw(z, x, M_signal, M_noise, BT_lag, use_BT)

% 
% [shat, thetahatfir] = p1_ncw_AR(z, x, M_signal, M_noise)
%
%   z           - Noisy sequence z(n) = s(n) + x(n)
%   x           - Noise samples x(n), when speaker is silent
%   M_signal    - AR model order to estimate z(n)
%   M_noise     - AR model order to estimate noise
%   BT_lag      - Blackman-Tuckey's method lag. Higher lag provides better
%                 spectral resolution but with higher variance
%   use_BT      - Binary. If '1', use BT non-parametric estimation. If '0', 
%                 use AR-M parametric estimation               
%
%   shat        - Estimate of s(n)
%   numnc,dennc - Non-causal filter
%
% Background noise cancellation with non-causal filter. Noise x(n) is
% estimated with an AR-M_noise model. Signal z(n) can be estimated either
% with a parametric method (AR-M_signal model), or non-parametric method
% (Blackman-Tuckey)
%
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get AR estimate for noise x(n)
    [Anoisehat, sigma2noisehat] = aryule(x, M_noise);
    [PhixxNum, PhixxDen] = filtspec(1, Anoisehat, sigma2noisehat);
    
    
    % -------- Parametric estimation of z(n) ---------
    % ---------------- AR-M model --------------------
    
    [Ahat, sigma2hat] = aryule(z, M_signal);
    [PhizzNum_AR, PhizzDen_AR] = filtspec(1, Ahat, sigma2hat);        
 
    
    % ------ Non-Parametric estimation of z(n) -------
    % ---------- Blackman-Tuckey's method ------------
    
    % Select window
    window = hamming(2*BT_lag-1);
    window = chebwin(2*BT_lag-1);
    window = blackman(2*BT_lag-1);

    % Compute BT spectrum estimate
    SigmaZzhat = xcovhat(z, z, BT_lag);
    PhizzNum_BT = [flip(SigmaZzhat(2:BT_lag)); SigmaZzhat];
    PhizzNum_BT = (PhizzNum_BT .* window)';
    PhizzDen_BT = [zeros(BT_lag-1,1); 1]';
   
    % -------------------------------------------------
           
    % Use Blackman-Tuckey estimation
    if use_BT
        PhizzNum = PhizzNum_BT;
        PhizzDen = PhizzDen_BT;
    else 
        PhizzNum = PhizzNum_AR;
        PhizzDen = PhizzDen_AR;
    end
    
    % Get cross-spectrum estimate Phizs
    % Since z(n) and x(n) are uncorrelated, Phizs = Phiss = Phizz - Phixx
    [PhissNum, PhissDen] = substract(PhizzNum, PhizzDen, PhixxNum, PhixxDen);
    PhizsNum = PhissNum;
    PhizsDen = PhissDen;
    
    % Compute Non-causal filter
    numnc = conv(PhizsNum, PhizzDen);
    dennc = conv(PhizsDen, PhizzNum);

    % Estimate s(n)
    shat = ncfilt(numnc, dennc, z);
    
    
    % --- PLOTS ---
    
    % Plot parametric vs non-parametric spectrum estimates
    w=linspace(0, pi, 512);
    [magp,~,wp]=dbode(PhizzNum_AR, PhizzDen_AR, 1, w);
    [magbt,~,wbt]=dbode(PhizzNum_BT, PhizzDen_BT, 1, w);
    figure
    plt = semilogy(wp, magp.^2, wbt, magbt.^2);
    set(plt, 'LineWidth', 1.5)
    legend('Parametric estimation (AR-30)', 'Non-parametric (Blackman-Tuckey)')
    title('Input z(n) PSD estimate')
    grid on
    set(gca,'FontSize', 14)
    set(legend,'Interpreter','latex')
    
    % Plot input, noise, output spectra and filter freq. response
    [Aouthat, sigma2outhat] = aryule(shat, M_signal);
    
    w=linspace(0, pi, 512);
    [magz,~,wz]=dbode(1,Ahat,1,w);
    [mags,~,ws]=dbode(1,Aouthat,1,w);
    [magx,~,wx]=dbode(1,Anoisehat,1,w);
    [magnc,~,wnc]=dbode(numnc,dennc,1,w);
    magout = magz.^2*sigma2hat .* magnc.^2;
    figure
    plt = semilogy(wz,magz.^2*sigma2hat, 'b', ...   % Input
                   ws, magout, 'r', ... % Output as Input * filter - should we use this, or AR estimate of output?
                   ... % ws, mags.^2*sigma2outhat, 'r', 
                   wx, magx.^2*sigma2noisehat, 'k:', ...
                   wnc, magnc.^2, '--');
    set(plt, 'LineWidth', 1.5)
    title('Non-causal Wiener fitler')
    legend('Input z(n) PSD estimate (AR-30)','Output $\hat{s}$(n)','Noise PSD estimate (AR-10)', 'Non-causal filter freq. response')
    set(legend,'Interpreter','latex')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
    grid on
    set(gca,'FontSize', 14)
%     
%     hold on
%     w=linspace(0, pi, 512);
%     [magz,~,wz]=dbode(1,Ahat,1,w);
%     [magnc,~,wnc]=dbode(numnc,dennc,1,w);
%     magout = magz.^2*sigma2hat .* magnc.^2;
%     %figure
%     plt = semilogy(wz, magout);
%     set(plt, 'LineWidth', 1.5)
%     legend('FIR filter (length 10)','FIR filter (length 30)', 'Causal filter','Non-causal filter')
%     set(legend,'Interpreter','latex')
%     xlabel('Frequency (rad/s)')
%     ylabel('Magnitude')
%     grid on
%     set(gca,'FontSize', 14)
    
    