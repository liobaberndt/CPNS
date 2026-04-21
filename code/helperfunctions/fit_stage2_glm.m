function mdl = fit_stage2_glm(X, y, stage2Type)
    y = y(:);
    if isempty(X) || size(X, 2) == 0
        switch stage2Type
            case 'linear'
                tbl = table(y, 'VariableNames', {'y'});
                mdl = fitglm(tbl, 'y ~ 1');
            case 'negbin'
                mdl = fit_stage2_negbin_matrix([], y);
            case 'logistic'
                tbl = table(y, 'VariableNames', {'y'});
                mdl = fitglm(tbl, 'y ~ 1', 'Distribution', 'Binomial');
            otherwise
                error('Unknown stage2Type: %s', stage2Type);
        end
        return;
    end

    switch stage2Type
        case 'linear'
            mdl = fitglm(X, y);
        case 'negbin'
            mdl = fit_stage2_negbin_matrix(X, y);
        case 'logistic'
            mdl = fitglm(X, y, 'Distribution', 'Binomial');
        otherwise
            error('Unknown stage2Type: %s', stage2Type);
    end
end

function mdl = fit_stage2_negbin_matrix(X, y)
    y = y(:);
    if any(~isfinite(y)) || any(y < 0) || any(abs(y - round(y)) > 0)
        error('fit_stage2_glm:NegBinY', 'Negative Binomial requires non-negative integer outcomes.');
    end

    if isempty(X)
        Z = ones(numel(y), 1);
    else
        Z = [ones(size(X, 1), 1), X];
    end

    beta0 = zeros(size(Z, 2), 1);
    ybar = mean(y);
    if ~isfinite(ybar) || ybar <= 0
        ybar = 0.1;
    end
    beta0(1) = log(ybar);
    logTheta0 = 0;
    p0 = [beta0; logTheta0];

    obj = @(p) negbin_nll(p, Z, y);
    opts = optimset('Display', 'off', 'MaxIter', 2e4, 'MaxFunEvals', 2e4, 'TolX', 1e-8, 'TolFun', 1e-8);
    phat = fminsearch(obj, p0, opts);

    beta = phat(1:end-1);
    logTheta = phat(end);
    theta = exp(logTheta);
    mu = exp(Z * beta);

    H = finite_diff_hessian(obj, phat);
    covP = pinv(H);
    seP = sqrt(max(diag(covP), 0));
    seBeta = seP(1:numel(beta));

    z = beta ./ max(seBeta, eps);
    pVals = 2 * (1 - normcdf(abs(z)));

    coefNames = cell(numel(beta), 1);
    coefNames{1} = '(Intercept)';
    for i = 2:numel(beta)
        coefNames{i} = sprintf('x%d', i - 1);
    end

    coefTable = table(beta, seBeta, z, pVals, 'VariableNames', {'Estimate', 'SE', 'tStat', 'pValue'}, ...
        'RowNames', coefNames);

    mdl = struct();
    mdl.model_class = 'custom_negbin_glm';
    mdl.link = 'log';
    mdl.dispersion_theta = theta;
    mdl.beta = beta;
    mdl.log_theta = logTheta;
    mdl.Coefficients = coefTable;
    mdl.Fitted = mu;
    mdl.NegativeLogLikelihood = obj(phat);
end

function nll = negbin_nll(p, Z, y)
    beta = p(1:end-1);
    logTheta = p(end);
    theta = exp(logTheta);
    eta = Z * beta;
    mu = exp(eta);
    mu = max(mu, eps);

    ll = gammaln(y + theta) - gammaln(theta) - gammaln(y + 1) ...
        + theta .* log(theta ./ (theta + mu)) ...
        + y .* log(mu ./ (theta + mu));
    nll = -sum(ll);
    if ~isfinite(nll)
        nll = realmax;
    end
end

function H = finite_diff_hessian(fun, x0)
    d = numel(x0);
    H = zeros(d, d);
    step = 1e-4 * (1 + abs(x0));
    f0 = fun(x0);
    for i = 1:d
        ei = zeros(d, 1);
        ei(i) = step(i);
        fpp = fun(x0 + ei);
        fmm = fun(x0 - ei);
        H(i, i) = (fpp - 2 * f0 + fmm) / (step(i)^2);
        for j = i+1:d
            ej = zeros(d, 1);
            ej(j) = step(j);
            f1 = fun(x0 + ei + ej);
            f2 = fun(x0 + ei - ej);
            f3 = fun(x0 - ei + ej);
            f4 = fun(x0 - ei - ej);
            Hij = (f1 - f2 - f3 + f4) / (4 * step(i) * step(j));
            H(i, j) = Hij;
            H(j, i) = Hij;
        end
    end
end
