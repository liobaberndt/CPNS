function [baseline_psd, source] = get_baseline_psd(state, receptor, baseline_data_path)
    f0 = fullfile(baseline_data_path, sprintf('%s_%s_Adjustment_0.00.mat', state, receptor));
    if ~isfile(f0)
        error('Baseline sweep missing (required): %s', f0);
    end
    data = load(f0, 'sim_psd');
    if ~isfield(data, 'sim_psd')
        error('Baseline file has no sim_psd: %s', f0);
    end
    baseline_psd = normalize_psd_vector(data.sim_psd, f0);
    source = 'adjustment_0';
end
