
function [xhat, num, den] = ncw(y, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen)
    % Apply non-causal Wiener filter (ncw) to input y to generate estimate
    % xhat - WORKS
    %
    % Model:
    %   x = H * y
    %   H = Phixy / Phiyy = PhixyNum·PhiyyDen / PhyxyDen·PhyyyNum
    %
    % output:
    %   xhat: estimate
    %   num: numerator of non-causal filter H
    %   den: denominator of non-causal filter H

    num = conv(PhixyNum, PhiyyDen);
    den = conv(PhixyDen, PhiyyNum);

    xhat = ncfilt(num, den, y);