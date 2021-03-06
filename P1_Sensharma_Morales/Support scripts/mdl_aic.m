function [] = mdl_aic(sig, M_lim)

sigma2array = zeros(1, M_lim);
aic_array = zeros(1, M_lim);
mdl_array = zeros(1, M_lim);

N = length(sig);

for M = 1:M_lim
    [~, sigma2hat] = aryule(sig, M);
    sigma2array(M) = log(sigma2hat);
    aic = log(sigma2hat) + (2 * M / N);
    mdl = (N * log(sigma2hat)) + (M * log(N));
    aic_array(M) = aic;
    mdl_array(M) = mdl;
end

M_count = linspace(1, M_lim, M_lim);
[mval_s2, midx_s2] = min(sigma2array)
[mval_aic, midx_aic] = min(aic_array)
[mval_mdl, midx_mdl] = min(mdl_array)

figure;
plot(M_count, sigma2array);
hold on;
plot(M_count, aic_array);
title('Best AR order estimation - AIC')
xlabel('AR Order (M)')
ylabel('Variance and AIC')
legend('Variance', 'AIC')

figure;
plot(M_count, mdl_array);
title('Best AR order estimation - MDL')
xlabel('AR Order (M)')
ylabel('MDL')
