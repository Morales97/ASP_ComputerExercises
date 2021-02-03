function [w, BT_spectrum] = BlackmanTuckey(y)
    % Blackman Tuckey spectral estimation method
    %
    % P(nu) = FT{r(k)·w(k)}
    % in freq. domain, performs a smoothing of the periodogram with W(nu)
    % 
    
    % Get ACF and multiply with window 
    M = 40;
    N = 1024;
    window = hamming(2*M+1);
    K = length(y);
    r_y = zeros(M+1, 1);
    for k = 0:M
        seq1 = y(1:K-k);
        seq2 = y(k+1:K);
        r_y(k+1) = sum(seq1 .* seq2);
    end
    r_y = r_y/K;
    r_y = [flip(r_y(2:M+1)); r_y];
    
    BT_spectrum = abs(fft(r_y .* window, N));
    w = linspace(0,pi,N/2);
    plot(w, BT_spectrum(1:N/2))