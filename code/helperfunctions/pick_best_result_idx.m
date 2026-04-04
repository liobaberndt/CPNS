function idx = pick_best_result_idx(entry)
    n = numel(entry.adjustments);
    if n == 0
        idx = 1;
        return;
    end
    if all(isnan(entry.mse_vs_sibling))
        [~, idx] = min(abs(entry.adjustments));
        return;
    end
    valid = isfinite(entry.mse_vs_sibling);
    if ~any(valid)
        [~, idx] = min(abs(entry.adjustments));
        return;
    end
    mmin = min(entry.mse_vs_sibling(valid));
    cand = find(valid & abs(entry.mse_vs_sibling - mmin) < 1e-12);
    if numel(cand) == 1
        idx = cand;
        return;
    end
    es = entry.effect_size(cand);
    if all(isnan(es))
        idx = cand(1);
        return;
    end
    [~, k] = max(es);
    idx = cand(k);
end
