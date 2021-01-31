
function [PhixyNum, PhixyDen, PhiyyNum, PhiyyDen] = spec_add(A, sigma2, Anoise, sigma2noise)
    % Find spectrum of y and cross-spectrum xy - WORKS
    %
    % Model:
    %   x is an AR process defined by A
    %   v is AR noise defined by Anoise
    %   y = x + v
    %
    % Input:
    %   A: AR polynomial
    %   sigma2: variance of white noise used to generate x
    %   Anoise: AR polynomial of noise
    %   sigma2noise: variance of white noise used to generate v
    %
    % Output:
    %   spectrum of x and cross-spectrum xy, divided into its num and dem

    [PhixxNum, PhixxDen] = filtspec(1, A, sigma2);
    [PhivvNum, PhivvDen] = filtspec(1, Anoise, sigma2noise);

    % x and v uncorrelated, y is addition of their sectra
    [PhiyyNum, PhiyyDen] = add(PhixxNum, PhixxDen, PhivvNum, PhivvDen);
    
    % x and v uncorrelated, cross-spectrum xy is spectrum of x
    PhixyNum = PhixxNum;
    PhixyDen = PhixxDen;