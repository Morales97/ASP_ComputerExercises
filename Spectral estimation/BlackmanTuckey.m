function [BT_spectrum] = BlackmanTuckey(y, window)
    % Blackman Tuckey spectral estimation method
    %
    % P(nu) = FT{r(k)·w(k)}
    % in freq. domain, performs a smoothing of the periodogram with W(nu)
    % 
    
    % Get ACF and multiply with window 
    K = length(y);
    N = length(window);
    M = ceil(N/2);
    r_y = zeros(M, 1);
    for k = 0:M-1
        seq1 = y(1:K-k);
        seq2 = y(k+1:K);
        r_y(k+1) = sum(seq1 .* seq2);
    end
    r_y = r_y/K;
    r_y = [flip(r_y(2:M)); r_y];
    
    BT_spectrum = 20 * log(abs(fft(r_y .* window, 1024)));