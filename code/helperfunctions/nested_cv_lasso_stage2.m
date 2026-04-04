function out = nested_cv_lasso_stage2(X, y, lassoDist, stage2Type, varargin)
% Nested CV LASSO + stage-2 GLM. Optional name-value: XCov (unpenalized columns; [] = all columns penalized).

    p = inputParser;
    addParameter(p, 'outerK', 5, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'innerK', 5, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'use1SE', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'minFoldFraction', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    addParameter(p, 'XCov', [], @(x) isnumeric(x) || isempty(x));
    parse(p, varargin{:});
    outerK = p.Results.outerK;
    innerK = p.Results.innerK;
    use1SE = logical(p.Results.use1SE);
    minFoldFraction = p.Results.minFoldFraction;
    XCov = p.Results.XCov;

    y = real(double(y(:)));
    y(~isfinite(y)) = NaN;
    [n, nPen] = size(X);
    if numel(y) ~= n && numel(y) == nPen
        X = X.';
        [n, nPen] = size(X);
    end
    if numel(y) ~= n
        error('nested_cv_lasso_stage2: X and y row count mismatch (size(X)=%s, numel(y)=%d).', mat2str(size(X)), numel(y));
    end

    partial = ~isempty(XCov);
    if partial
        if size(XCov, 1) ~= n
            error('nested_cv_lasso_stage2:XCovRows', 'XCov must have same number of rows as X (%d).', n);
        end
        nCov = size(XCov, 2);
        ok = ~any(isnan(X), 2) & ~any(isnan(XCov), 2) & ~isnan(y) & isfinite(y);
    else
        nCov = 0;
        ok = ~any(isnan(X), 2) & ~isnan(y) & isfinite(y);
    end
    if ~all(ok)
        X = X(ok, :);
        y = y(ok);
        n = size(X, 1);
        if partial
            XCov = XCov(ok, :);
        end
    end

    if n < outerK + 2
        error('nested_cv_lasso_stage2: need more observations (%d) than outer folds (%d).', n, outerK);
    end

    cvo = cvpartition(n, 'KFold', outerK);
    yhat = nan(n, 1);
    selection_counts = zeros(1, nPen);
    cv_mae = zeros(outerK, 1);

    for k = 1:outerK
        tr = training(cvo, k);
        te = test(cvo, k);
        Xtr = X(tr, :);
        ytr = y(tr);
        Xte = X(te, :);
        if partial
            Xtr_cov = XCov(tr, :);
            Xte_cov = XCov(te, :);
        end

        if sum(tr) < innerK + 1
            error('nested_cv_lasso_stage2: training fold too small for inner CV.');
        end

        idxL = [];
        B = [];
        FitInfo = [];

        if strcmp(lassoDist, 'gaussian')
            [B, FitInfo] = lasso(Xtr, ytr, 'CV', innerK);
            if use1SE && isfield(FitInfo, 'Index1SE')
                idxL = FitInfo.Index1SE;
            else
                idxL = FitInfo.IndexMinMSE;
            end
            lam = FitInfo.Lambda(idxL);
            [B, FitInfo] = lasso(Xtr, ytr, 'Lambda', lam);
            active = abs(B) > 1e-10;
            if partial
                Xtr_s = [Xtr_cov, Xtr(:, active)];
                Xte_s = [Xte_cov, Xte(:, active)];
            else
                Xtr_s = Xtr(:, active);
                Xte_s = Xte(:, active);
            end
        else
            if strcmp(lassoDist, 'binomial')
                [B, FitInfo] = lassoglm(Xtr, ytr, 'binomial', 'CV', innerK, 'Alpha', 1);
            else
                [B, FitInfo] = lassoglm(Xtr, ytr, lassoDist, 'CV', innerK, 'Alpha', 1);
            end
            if use1SE && isfield(FitInfo, 'Index1SE')
                idxL = FitInfo.Index1SE;
            elseif isfield(FitInfo, 'IndexMinDeviance')
                idxL = FitInfo.IndexMinDeviance;
            else
                idxL = FitInfo.IndexMinMSE;
            end
            lam = FitInfo.Lambda(idxL);
            if strcmp(lassoDist, 'binomial')
                [B, FitInfo] = lassoglm(Xtr, ytr, 'binomial', 'Lambda', lam, 'Alpha', 1);
            else
                [B, FitInfo] = lassoglm(Xtr, ytr, lassoDist, 'Lambda', lam, 'Alpha', 1);
            end
            Bcol = B(:, end);
            beta = nested_cv_lassoglm_beta_from_B(Bcol, nPen);
            active = abs(beta) > 1e-10;
            active = active(:).';
            if partial
                Xtr_s = [Xtr_cov, Xtr(:, active)];
                Xte_s = [Xte_cov, Xte(:, active)];
            else
                Xtr_s = Xtr(:, active);
                Xte_s = Xte(:, active);
            end
        end

        selection_counts = selection_counts + active(:).';

        try
            mdl = fit_stage2_glm(Xtr_s, ytr, stage2Type);
            yhat_te = predict(mdl, Xte_s);
            yhat_te = yhat_te(:);
        catch
            yhat_te = nan(numel(y(te)), 1);
        end

        yhat(te) = yhat_te;
        yte = y(te);
        yte = yte(:);
        cv_mae(k) = mean(abs(yte - yhat_te), 'omitnan');
    end

    valid = ~isnan(yhat) & ~isnan(y);
    if numel(y(valid)) >= 3
        cv_r = corr(y(valid), yhat(valid), 'rows', 'complete');
    else
        cv_r = NaN;
    end

    minVotes = max(1, ceil(outerK * minFoldFraction));
    active_pen = selection_counts >= minVotes;
    active_pen = active_pen(:).';

    if partial && ~any(active_pen)
        warning('nested_cv_lasso_stage2:NoStablePen', ...
            ['No DCM-related predictor reached >= %d/%d outer folds; ', ...
            'stage-2 uses forced covariates only.'], minVotes, outerK);
    elseif ~partial && ~any(active_pen)
        warning('nested_cv_lasso_stage2:NoStableSelection', ...
            'No predictor reached >= %d/%d outer folds; intercept-only stage-2.', minVotes, outerK);
    end

    if partial
        Xs = [XCov, X(:, active_pen)];
    else
        Xs = X(:, active_pen);
    end

    try
        final_mdl = fit_stage2_glm(Xs, y, stage2Type);
    catch
        final_mdl = [];
    end

    if partial
        pred_idx_final = [ (1:nCov)'; nCov + find(active_pen)' ];
    else
        pred_idx_final = find(active_pen);
    end

    out = struct();
    out.yhat_cv = yhat;
    out.cv_r = cv_r;
    out.cv_mae_mean = mean(cv_mae);
    out.selection_counts = selection_counts;
    out.active_pen = active_pen;
    out.active_final = active_pen;
    out.pred_idx_final = pred_idx_final(:);
    out.final_mdl = final_mdl;
    out.lassoDist = lassoDist;
    out.stage2Type = stage2Type;
    out.min_fold_fraction = minFoldFraction;
    out.min_votes = minVotes;
    out.n_cov = nCov;
end

function beta = nested_cv_lassoglm_beta_from_B(Bcol, nPred)
    k = numel(Bcol);
    if k == nPred
        beta = Bcol;
    elseif k == nPred + 1
        beta = Bcol(2:end);
    else
        error('nested_cv_lasso_stage2: lassoglm coefficient length %d does not match nPred=%d or nPred+1.', k, nPred);
    end
end
