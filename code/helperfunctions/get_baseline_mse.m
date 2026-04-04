function baseline_mse = get_baseline_mse(state, receptor, baseline_data_path, sibling_psd)
    [baseline_psd, ~] = get_baseline_psd(state, receptor, baseline_data_path);
    [base_cmp, sibling_cmp] = align_psd_vectors_or_error(baseline_psd, sibling_psd, ...
        sprintf('%s_%s_Adjustment_0.00.mat', state, receptor), 'baseline', 'sibling');
    diff_baseline = base_cmp - sibling_cmp;
    baseline_mse = mean(diff_baseline.^2);
end
