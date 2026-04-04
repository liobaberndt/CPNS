function mdl = fit_stage2_glm(X, y, stage2Type)
    y = y(:);
    n = numel(y);
    if isempty(X) || size(X, 2) == 0
        tbl = table(y, 'VariableNames', {'y'});
        switch stage2Type
            case 'linear'
                mdl = fitglm(tbl, 'y ~ 1');
            case 'negbin'
                mdl = fit_stage2_negbin_intercept(tbl);
            case 'logistic'
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

function mdl = fit_stage2_negbin_intercept(tbl)
    try
        mdl = fitglm(tbl, 'y ~ 1', 'Distribution', 'NegativeBinomial');
    catch
        try
            mdl = fitglm(tbl, 'y ~ 1', 'Distribution', 'Poisson');
        catch
            mdl = fitglm(tbl, 'y ~ 1');
        end
    end
end

function mdl = fit_stage2_negbin_matrix(X, y)
    try
        mdl = fitglm(X, y, 'Distribution', 'NegativeBinomial');
    catch
        try
            mdl = fitglm(X, y, 'Distribution', 'Poisson');
        catch
            mdl = fitglm(X, y);
        end
    end
end
