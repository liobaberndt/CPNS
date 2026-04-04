% Summarise neurotransmitter receptor sweep simulations against sibling reference PSDs.
%
% Input filenames: <Stage>_<Receptor>_Adjustment_<pct>.mat
%
% Workflow summary:
%   1) Configure paths, sleep stages, receptors, permutation settings
%   2) Load reference sibling PSDs and verify baseline (0%) sweep files exist
%   3) Scan folder, parse filenames, accumulate MSE/effect size per adjustment
%   4) Deduplicate adjustments, build summary (best adjustment, permutation p, BH-FDR q)
%   5) Save MAT/CSV

%==========================================================================
% SECTION 1: PATHS AND CONFIGURATION
%==========================================================================
clearvars

scriptPath = fileparts(mfilename('fullpath'));
cd(scriptPath);
addpath(scriptPath);
addpath(fullfile(scriptPath, 'helperfunctions'));

simulation_data_path = '/Users/liobaberndt/Library/CloudStorage/Dropbox/UoE/sleepdetectives/dcm/eLifeNick/DCM_v01/Simulations_4/';
baseline_data_path   = simulation_data_path;
reference_data_path  = simulation_data_path;
output_path          = simulation_data_path;

states = {'N2'};
receptor_types = {'AMPA', 'GABAA', 'GABAB', 'NMDA'};
num_permutations = 1000;
rng_seed = 1;
rng(rng_seed, 'twister');

fprintf('Simulation data path: %s\n', simulation_data_path);
fprintf('Reference data path: %s\n', reference_data_path);
fprintf('Output path: %s\n', output_path);

%==========================================================================
% SECTION 2: REFERENCE PSDS AND BASELINE SWEEP FILES
%==========================================================================
reference = struct();
for s = 1:length(states)
    state = states{s};
    reference.(state) = struct('sibling', []);
    sibling_file = fullfile(reference_data_path, sprintf('avg_sibling_%s.mat', state));
    if ~isfile(sibling_file)
        error('Required sibling reference missing: %s', sibling_file);
    end
    reference.(state).sibling = load_reference_psd(sibling_file, {'avg_sibling_psd', 'sim_psd'});
end

for s = 1:length(states)
    state = states{s};
    for r = 1:length(receptor_types)
        receptor = receptor_types{r};
        b0 = fullfile(baseline_data_path, sprintf('%s_%s_Adjustment_0.00.mat', state, receptor));
        if ~isfile(b0)
            error('Required baseline sweep missing: %s', b0);
        end
    end
end

%==========================================================================
% SECTION 3: LOAD SWEEP FILES AND BUILD RESULTS
%==========================================================================
files = dir(fullfile(simulation_data_path, '*_Adjustment_*.mat'));
if isempty(files)
    error('No simulation sweep files found in %s', simulation_data_path);
end
fprintf('Candidate sweep files: %d in %s\n', length(files), simulation_data_path);

%----------------------------------------------------------------------
% 3.1 Initialise per-stage / per-receptor result buckets
%----------------------------------------------------------------------
results = struct();
for s = 1:length(states)
    state = states{s};
    results.(state) = struct();
    for r = 1:length(receptor_types)
        receptor = receptor_types{r};
        results.(state).(receptor) = struct( ...
            'adjustments', [], ...
            'mse_vs_sibling', [], ...
            'effect_size', [], ...
            'sim_psd_rows', {{}}, ...
            'source_file', {{}});
    end
end

%----------------------------------------------------------------------
% 3.2 Parse each file, MSE vs sibling, effect size vs 0% baseline
%----------------------------------------------------------------------
for i = 1:length(files)
    file_name = files(i).name;
    [ok, state, receptor, adjustment] = parse_sim_filename(file_name, states, receptor_types);
    if ~ok
        continue;
    end

    loaded = load(fullfile(files(i).folder, file_name), 'sim_psd');
    if ~isfield(loaded, 'sim_psd')
        error('File has no variable sim_psd: %s', file_name);
    end
    sim_psd = normalize_psd_vector(loaded.sim_psd, file_name);

    sibling_psd = reference.(state).sibling;
    [sim_cmp, sibling_cmp] = align_psd_vectors_or_error(sim_psd, sibling_psd, file_name, 'sim', 'sibling');
    diff_sim = sim_cmp - sibling_cmp;
    metric = mean(diff_sim.^2);

    baseline_mse = get_baseline_mse(state, receptor, baseline_data_path, sibling_psd);
    es = NaN;
    if ~isnan(baseline_mse) && baseline_mse ~= 0
        es = (baseline_mse - metric) / baseline_mse;
    end

    bucket = results.(state).(receptor);
    bucket.adjustments(end+1, 1) = adjustment;
    bucket.mse_vs_sibling(end+1, 1) = metric;
    bucket.effect_size(end+1, 1) = es;
    bucket.sim_psd_rows{end+1, 1} = sim_psd;
    bucket.source_file{end+1, 1} = file_name;
    results.(state).(receptor) = bucket;
end

for s = 1:length(states)
    state = states{s};
    for r = 1:length(receptor_types)
        receptor = receptor_types{r};
        if isempty(results.(state).(receptor).adjustments)
            error('No sweep files matched %s / %s (check ''states'' vs filename prefixes in %s).', ...
                state, receptor, simulation_data_path);
        end
    end
end

%==========================================================================
% SECTION 4: DEDUPLICATE AND REPORT COVERAGE
%==========================================================================
for s = 1:length(states)
    state = states{s};
    fprintf('Coverage for %s:\n', state);
    for r = 1:length(receptor_types)
        receptor = receptor_types{r};
        entry = results.(state).(receptor);
        entry = deduplicate_adjustments_or_error(entry, state, receptor);
        results.(state).(receptor) = entry;
        n_unique = numel(unique(round(entry.adjustments, 6)));
        fprintf('  %s: %d unique adjustments\n', receptor, n_unique);
    end
end

%==========================================================================
% SECTION 5: SUMMARY TABLE (BEST ADJUSTMENT, PERMUTATION P, BH-FDR)
%==========================================================================
summary_rows = {};
summary_header = {'Sleep Stage', 'Receptor', 'N_Adjustments', 'Best Adjustment (%)', 'Best MSE', 'Effect Size', 'P-value', 'Q-value (BH-FDR)', 'Baseline Source', 'Best_Source_File', 'Sim_SumFirst10'};
perm_tests_executed = 0;
perm_tests_skipped_single_adjustment = 0;

for s = 1:length(states)
    state = states{s};
    for r = 1:length(receptor_types)
        receptor = receptor_types{r};
        entry = results.(state).(receptor);
        if isempty(entry.adjustments)
            continue;
        end
        unique_adjustments = unique(round(entry.adjustments, 6));
        n_adjustments = numel(unique_adjustments);
        if n_adjustments <= 1
            fprintf('WARNING: %s/%s has %d unique adjustment value (baseline-only or missing sweep).\n', ...
                state, receptor, n_adjustments);
        end

        idx = pick_best_result_idx(entry);
        if all(isnan(entry.mse_vs_sibling))
            best_mse = NaN;
        else
            best_mse = entry.mse_vs_sibling(idx);
        end

        best_adj = entry.adjustments(idx);
        best_es = entry.effect_size(idx);
        best_p = NaN;
        sibling_psd = reference.(state).sibling;
        [baseline_psd, baseline_source] = get_baseline_psd(state, receptor, baseline_data_path);
        if n_adjustments <= 1
            perm_tests_skipped_single_adjustment = perm_tests_skipped_single_adjustment + 1;
        else
            best_p = compute_permutation_pvalue(entry.sim_psd_rows{idx}, baseline_psd, sibling_psd, num_permutations);
            perm_tests_executed = perm_tests_executed + 1;
        end
        best_file = entry.source_file{idx};
        best_sim = entry.sim_psd_rows{idx};
        sim_tag = sim_digest(best_sim);
        summary_rows(end+1, :) = {state, receptor, n_adjustments, best_adj, best_mse, best_es, best_p, NaN, baseline_source, best_file, sim_tag};
    end
end

summary_rows = add_fdr_qvalues(summary_rows);

%==========================================================================
% SECTION 6: SAVE OUTPUTS
%==========================================================================
if ~isfolder(output_path)
    mkdir(output_path);
    fprintf('Created output directory: %s\n', output_path);
end

summary_mat = fullfile(output_path, 'sim_analysis_summary.mat');
save(summary_mat, 'results', 'summary_rows', 'summary_header');

summary_csv = fullfile(output_path, 'sim_analysis_summary.csv');
write_summary_csv(summary_csv, summary_header, summary_rows);

fprintf('Saved summary MAT: %s\n', summary_mat);
fprintf('Saved summary CSV: %s\n', summary_csv);
fprintf('Permutation tests executed: %d (x %d each = %d total permutations)\n', ...
    perm_tests_executed, num_permutations, perm_tests_executed * num_permutations);
fprintf('Permutation tests skipped (single adjustment only): %d\n', perm_tests_skipped_single_adjustment);
fprintf('Random seed used for permutation test: %d\n', rng_seed);