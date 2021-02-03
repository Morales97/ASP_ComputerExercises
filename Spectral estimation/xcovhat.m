function SigmaYxhat = xcovhat(x,y,N)
    % Find estimate of cross-covariance vector Nx1, given samples from x
    % and y
    %
    % Input
    %   x: process
    %   y: observations
    %   N: length of matrix NxN
    % 
    % Ouput
    %   SigmaYxhat: estimate of cross-covariance vector 

    r_yx = zeros(N,1);
    M = min(length(x),length(y));
    for k = 0:N-1
        seq1 = y(1:M-k);
        seq2 = x(k+1:M);
        r_yx(k+1) = sum(seq1 .* seq2);
    end
    SigmaYxhat = r_yx/M;
    
    