function SigmaYxhat = xcovhat(x,y,N)
%
% SigmaYxhat = xcovhat(x, y, N)
%
%	y, x			- Data sequences
%	N			- Size of covariance matrix (Nx1)
%
%  xcovhat: Estimates SigmaYx=E[Y(n)x(n)]
%
%		where 
%
%	   	Y(n)=[y(n) y(n-1) ... y(n-N+1)]^{T}
%
%     
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    r_yx = zeros(N,1);
    M = min(length(x),length(y));
    for k = 0:N-1
        seq1 = y(1:M-k);
        seq2 = x(k+1:M);
        r_yx(k+1) = sum(seq1 .* seq2);
    end
    SigmaYxhat = r_yx/M;
    
    