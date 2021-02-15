function SigmaYYhat = covhat(y, N)

%
% SigmaYYhat = covhat(y,N)
%
%	y			- Data sequence
%	N			- Size of covariance matrix
%
%  covhat: Estimates SigmaYY=E[Y(n)Y^{T}(n)]
%
%		where 
%
%	   	Y(n)=[y(n) y(n-1) ... y(n-N+1)]^{T}
%
%     
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1. Find BIASED acf estimator
    r_y = zeros(1, N);
    M = length(y);
    for k = 0:N-1
        seq1 = y(1:M-k);
        seq2 = y(k+1:M);
        r_y(k+1) = sum(seq1 .* seq2);
    end
    r_y = r_y/M;
    
    % 2. Generate covariance matrix from acf estimation (toeplitz)
    SigmaYYhat = toeplitz(r_y);