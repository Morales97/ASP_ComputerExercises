function [shat, numc, denc] = p1_cw(z, x, M_signal, M_noise)    

% 
% [shat, thetahatfir] = p1_cw(z, x, M_signal, M_noise)
%
%   z           - Noisy sequence z(n) = s(n) + x(n)
%   x           - Noise samples x(n), when speaker is silent
%   M_signal    - AR model order to estimate z(n)
%   M_noise     - AR model order to estimate noise
%
%   shat        - Estimate of s(n)
%   numc, denc  - Causal filter
%
% Background noise cancellation with causal filter.
%
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get AR estimate for noise x(n)
    [Anoisehat, sigma2noisehat] = aryule(x, M_noise);
    [PhixxNum, PhixxDen] = filtspec(1, Anoisehat, sigma2noisehat);
    
    % Get AR estimate for z(n)
    [Ahat, sigma2hat] = aryule(z, M_signal);
    [PhizzNum, PhizzDen] = filtspec(1, Ahat, sigma2hat);  
    
    % Get cross-spectrum estimate Phizs
    % Since z(n) and x(n) are uncorrelated, Phizs = Phiss = Phizz - Phixx
    [PhissNum, PhissDen] = substract(PhizzNum, PhizzDen, PhixxNum, PhixxDen);
    PhizsNum = PhissNum;
    PhizsDen = PhissDen;
    
    % Compute causal filter and estimate of s(n)
    % Select delay m=0, since we are filtering
    [shat, numc, denc] = cw(z, PhizsNum, PhizsDen, PhizzNum, PhizzDen, 0);
    
    
    % --- PLOTS ---
    [Aouthat, sigma2outhat] = aryule(shat, M_signal);
    
    figure
    w=linspace(0, pi, 512);
    [magz,~,wz]=dbode(1,Ahat,1,w);
    [mags,~,ws]=dbode(1,Aouthat,1,w);
    [magx,~,wx]=dbode(1,Anoisehat,1,w);
    [magc,~,wc]=dbode(numc,denc,1,w);
    magout = magz.^2*sigma2hat .* magc.^2;
    plt = semilogy(wz, magz.^2*sigma2hat, 'b', ...
                   ws, magout, 'r', ...
                   wx, magx.^2*sigma2noisehat, ':k', ...
                   wc, magc.^2, '--');
    set(plt, 'LineWidth', 1.5)
    title('Causal Wiener filtering')
    legend('Input z(n) PSD estimate (AR-30)','Output $\hat{s}$(n)','Noise PSD estimate (AR-10)', 'Causal filter freq. response')
    set(legend,'Interpreter','latex')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
    grid on
    set(gca,'FontSize', 14)
    
%     % Comparison
%     hold on
%     w=linspace(0, pi, 512);
%     [magz,~,wz]=dbode(1,Ahat,1,w);
%     [magc,~,wc]=dbode(numc,denc,1,w);
%     magout = magz.^2*sigma2hat .* magc.^2;
%     %figure
%     plt = semilogy(wz, magout);
%     set(plt, 'LineWidth', 1.5)
%     legend('Causal filter')
%     set(legend,'Interpreter','latex')
%     xlabel('Frequency (rad/s)')
%     ylabel('Magnitude')
%     grid on
%     set(gca,'FontSize', 14)
    