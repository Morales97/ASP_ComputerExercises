function SigmaYYhat = covhat(y, N)
    % Find estimate of covariance matrix NxN, given sequence samples
    %
    % Steps
    % 1. Find BIASED acf estimator
    % 2. Generate covariance matrix from acf estimation (toeplitz)
    %
    % Input
    %   y: sequence
    %   N: length of matrix NxN
    % 
    % Ouput
    %   SigmaYYhat: estimate of covariance matrix 
    %   ATTENTION - r_y(0) will be unbiased, but for all k>0 will be
    %   biased. The difference should be very small for large M
    
    % Get biased acf estimator of length N
    r_y = zeros(1, N);
    M = length(y);
    for k = 0:N-1
        seq1 = y(1:M-k);
        seq2 = y(k+1:M);
        r_y(k+1) = sum(seq1 .* seq2);
    end
    r_y = r_y/M;
    
    % Generate covariance matrix as toeplitz
    SigmaYYhat = toeplitz(r_y);