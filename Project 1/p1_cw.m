function [shat, numc, denc] = p1_cw(z, x, N, ar_order)    
    % Project 1 - Causal Wiener filter
    %
    % Input
    %   z: commentator + vuvuzelas
    %   x: vuvuzelas (z(1:8000))
    %   N: number of autocorrelation samples used to get spectrum of z
    %   ar_order: order of the AR-m model that tries to fit into x
    
    % Get noise spectrum
    [Ahat, sigma2hat] = ar_id(x, ar_order);
    SigmaXxhat = xcovhat(x,x,N);
    [PhixxNum, PhixxDen] = filtspec(1, Ahat, sigma2hat);
    
    % Get Z spectrum
    SigmaZzhat = xcovhat(z,z,N);
    PhizzNum = [flip(SigmaZzhat(2:N)); SigmaZzhat]';
    PhizzDen = [zeros(N-1,1); 1]';
    
    PhizzNum_causal = SigmaZzhat;
    PhizzDen_causal = 1;
    
    PhizzNum_anticausal = flip(SigmaZzhat(2:N));
    PhizzDen_anticausal = 1;
    
    % Estimate original signal S spectrum
    SigmaSshat = SigmaZzhat - SigmaXxhat;
    PhissNum_causal = SigmaSshat;
    PhissDen_causal = 1;
    
    % Get causal part of Phiss / Phizz_anticausal
    num = PhissNum_causal;
    den = PhizzNum_anticausal;
        % Convert into num/den * delay. I think it's not necessary bc delay
        % will always be zero
    [num, den, dly] = delay(num,den);
        % Factorize denominator into causal / anticausal
    [den_causal, den_anticausal] = fact(den);
    
    % Filter
    numc = num;
    denc = conv(PhizzNum_causal, den_causal)';
    
    shat = filter(numc, denc, z);
    
    
    % Plot
    [wbt, BT_spectrum_z] = BlackmanTuckey(z);
    [wbt, BT_spectrum_shat] = BlackmanTuckey(shat);
    
    w=linspace(0,pi, 512);
    [magv,phasev,wv]=dbode(1,Ahat',1,w);
    [magc,~,wc]=dbode(numc,denc,1,w);
    plt = semilogy(wbt, BT_spectrum_z(1:512).^2, wbt, BT_spectrum_shat(1:512).^2, wv, magv.^2*sigma2hat, wc, magc.^2, '--');
    set(plt, 'LineWidth', 1.5)
    title('Spectra')
    legend('Input z (BT)', 'Output shat (BT)','Noise', 'Causal freq response')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
