function p_value = compute_permutation_pvalue(sim_psd, baseline_psd, sibling_psd, num_permutations)
    p_value = NaN;
    if num_permutations < 1
        return;
    end
    ns = numel(sim_psd);
    nb = numel(baseline_psd);
    nk = numel(sibling_psd);
    n = min([ns, nb, nk]);
    if n == 0
        return;
    end
    if ~(ns == nb && nb == nk)
        error('Permutation test requires equal-length PSDs (sim=%d, baseline=%d, sibling=%d).', ns, nb, nk);
    end
    sim_c = sim_psd(1:n);
    base_c = baseline_psd(1:n);
    sib_c = sibling_psd(1:n);

    ds = sim_c - sib_c;
    db = base_c - sib_c;
    mse_sim = mean(ds.^2);
    mse_base = mean(db.^2);
    if mse_base == 0
        return;
    end

    observed_es = (mse_base - mse_sim) / mse_base;
    perm_es = zeros(num_permutations, 1);
    for perm = 1:num_permutations
        swap = rand(n, 1) > 0.5;
        sim_p = ds(:);
        base_p = db(:);
        sim_p(swap) = db(swap);
        base_p(swap) = ds(swap);
        mse_sim_p = mean(sim_p.^2);
        mse_base_p = mean(base_p.^2);
        if mse_base_p == 0
            perm_es(perm) = 0;
        else
            perm_es(perm) = (mse_base_p - mse_sim_p) / mse_base_p;
        end
    end
    p_value = (sum(perm_es >= observed_es) + 1) / (num_permutations + 1);
end
