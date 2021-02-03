function [shat, numnc, dennc] = p1_ncw(z, x, N, ar_order)
    % Project 1 - Non-causal Wiener filter
    %
    % Input
    %   z: commentator + vuvuzelas
    %   x: vuvuzelas (z(1:8000))
    %   N: number of autocorrelation samples used to get spectrum of z
    %   ar_order: order of the AR-m model that tries to fit into x

    % Get noise spectrum
    [Ahat, sigma2hat] = ar_id(x, ar_order);
    [PhixxNum, PhixxDen] = filtspec(1, Ahat, sigma2hat);
    
    w=linspace(0,pi, 512);
    [mag1,~,wv]=dbode(1,Ahat',1,w);
    [mag2,~,wv]=dbode(PhixxNum,PhixxDen,1,w);
    
    semilogy(wv, mag1.^2*sigma2hat)
    figure
    semilogy(wv, mag2)
    figure
    % WARNING - for some reason, this plots are not the same (they should
    % be). Same curve, but differ by some factor  
    
    % Get signal spectrum
    SigmaZzhat = xcovhat(z,z,N);
    PhizzNum = [flip(SigmaZzhat(2:N)); SigmaZzhat]';
    PhizzDen = [zeros(N-1,1); 1]';
    [mag3,~,wv]=dbode(PhizzNum,PhizzDen,1,w);
    semilogy(wv, mag3.^2)
    
    % Filter
    numnc = conv(PhixxNum, PhizzDen);
    dennc = conv(PhixxDen, PhizzNum);

    shat = ncfilt(numnc, dennc, z);
    
    % plot
    [wbt, BT_spectrum_z] = BlackmanTuckey(z);
    [wbt, BT_spectrum_shat] = BlackmanTuckey(shat);
    
    w=linspace(0,pi, 512);
    [magv,phasev,wv]=dbode(1,Ahat',1,w);
    [magnc,~,wnc]=dbode(numnc,dennc,1,w);
    plt = semilogy(wbt, BT_spectrum_z(1:512).^2, wbt, BT_spectrum_shat(1:512).^2, wv, magv.^2*sigma2hat, wnc, magnc.^2, '--');
    set(plt, 'LineWidth', 1.5)
    title('Spectra')
    legend('Input z (BT)', 'Output shat (BT)','Noise', 'Non-causal freq response')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')