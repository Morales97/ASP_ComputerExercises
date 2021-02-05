function [shat, thetahatfir] = p1_firw(z, x, N, ar_order)
    % Background noise reduction with FIR Wiener filter
    %
    % Model
    %   z(n) = s(n) + x(n)
    %   z is the measurement of signal s + noise x
    %
    % Input
    %   z: noisy measurments
    %   x: vuvuzelas (z(1:8000))
    %   N: length of FIR filter
    %
    % Output
    %   shat: estimation of signal s
    %   theta: FIR parameters

    % 1. Get SigmaZZ
    SigmaZZhat = covhat(z, N);
    
    % 2. Get SigmaZs
    % SigmaZs = SigmaSs - SigmaXs = SigmaSs
    % SigmaZz = SigmaSs + SigmaXx 
    SigmaXxhat = xcovhat(x,x,N);
    % HEURISTICS - check if using samples with speech delivers better
    % results
    z_speaking = z(12000:30000);
    SigmaZzhat = xcovhat(z_speaking,z_speaking,N);
    %SigmaZzhat = xcovhat(z,z,N);
    SigmaSshat = SigmaZzhat - SigmaXxhat;
    SigmaZshat = SigmaSshat;
    
    % 3. Obtain filter coefficients
    thetahatfir = SigmaZZhat\SigmaZshat;
    
    % 4. Obtain filtered signal
    shat = conv(z, thetahatfir);

    % --- PLOT ---
    % Get spectrum using Blackman Tuckey's method
    [wbt, BT_spectrum_z] = BlackmanTuckey(z, 20);
    [wbt, BT_spectrum_shat] = BlackmanTuckey(shat, 20);
    % Get AR
    [Ahat, sigma2hat] = ar_id(x, ar_order);

    w=linspace(0,pi, 512);
    [magv,phasev,wv]=dbode(1,Ahat',1,w);
    [magfir,phasefir,wfir]=dbode(thetahatfir',1,1,w);
    plt = semilogy(wbt, BT_spectrum_z(1:512).^2, wbt, BT_spectrum_shat(1:512).^2, wv, magv.^2*sigma2hat, wfir, magfir.^2, '--');
    set(plt, 'LineWidth', 1.5)
    title('Spectra')
    legend('Input z (BT)', 'Output shat (BT)','Noise', 'FIR freq response')
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
