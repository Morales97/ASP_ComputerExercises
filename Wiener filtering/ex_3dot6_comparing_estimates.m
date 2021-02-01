% Load wrefdata before executing

% Compute FIR, non-Causal and Causal filters and estimates
[xhatnc, xhatc, xhatfir, numnc, dennc, numc, denc, thetahatfir] = est_add(x, v, N, A, sigma2, Anoise, sigma2noise, SigmaYx, SigmaYY)

% Compare spectra
figure()
spec_comp(A,sigma2,Anoise,sigma2noise,numnc,dennc,numc,denc,thetahatfir)