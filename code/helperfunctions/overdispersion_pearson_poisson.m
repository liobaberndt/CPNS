function [dispersion_ratio, p_gof] = overdispersion_pearson_poisson(y, mu)
% Pearson dispersion for Poisson fit and rough goodness-of-fit p-value.
% dispersion_ratio >> 1 suggests overdispersion vs Poisson (NB may be better).
    y = y(:);
    mu = max(mu(:), eps);
    pearson = sum(((y - mu) ./ sqrt(mu)).^2);
    df = max(numel(y) - 1, 1);
    dispersion_ratio = pearson / df;
    p_gof = 1 - chi2cdf(pearson, df);
end
