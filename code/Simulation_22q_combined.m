% Run neurotransmitter simulation sweeps on 22q sleep-stage DCMs.
%
% This script loads stage-specific DCM fits, runs parameter sweeps for
% AMPA/GABA-A/GABA-B/NMDA connection groups, and saves the resulting
% simulated patient spectra. Each sweep uses that subject's fitted DCM.Ep;
% saved .mat files contain the group-mean PSD across successful runs.
%
% Workflow summary:
%   1) Configure paths and parameter sets
%   2) Validate available stage files
%   3) Run neurotransmitter-wise parameter sweeps
%   4) Save simulation outputs

%==========================================================================
% SECTION 1: PATHS AND SETUP
%==========================================================================
clearvars

scriptPath = fileparts(mfilename('fullpath'));
cd(scriptPath);
addpath(scriptPath);
addpath(fullfile(scriptPath, 'helperfunctions'));

base_path            = '/Users/liobaberndt/Dropbox/UoE';
DCM_version          = '/Users/liobaberndt/Library/CloudStorage/Dropbox/UoE/sleepdetectives/dcm/eLifeNick/DCM_v01/';
atcm_dir             = fullfile(base_path, 'DCM_TCM', 'atcm');
simulation_save_path = fullfile(DCM_version, 'Simulations_test');

if ~isfolder(simulation_save_path)
    mkdir(simulation_save_path);
end

% Initialise SPM defaults.
spm('defaults', 'eeg');

%==========================================================================
% SECTION 2: NEUROTRANSMITTER PARAMETER DEFINITIONS
%==========================================================================
% Note: fixed-prior connections from RunTCM are excluded from sweep lists.

% AMPA Parameters
AMPA_params = {
    struct('name', 'AMPA_ss_sp', 'matrix', 'H',  'i', 2, 'j', 1),
    struct('name', 'AMPA_ss_rl', 'matrix', 'H',  'i', 8, 'j', 1),
    struct('name', 'AMPA_sp_si', 'matrix', 'H',  'i', 3, 'j', 2),
    struct('name', 'AMPA_sp_dp', 'matrix', 'H',  'i', 4, 'j', 2),
    struct('name', 'AMPA_si_sp', 'matrix', 'H',  'i', 2, 'j', 3),
    struct('name', 'AMPA_dp_tp', 'matrix', 'H',  'i', 6, 'j', 4),
    struct('name', 'AMPA_rl_ss', 'matrix', 'H',  'i', 1, 'j', 8)
};

% GABA-A Parameters
GABAA_params = {
    struct('name', 'GABAA_ss_ss', 'matrix', 'H', 'i', 1, 'j', 1),
    struct('name', 'GABAA_ss_sp', 'matrix', 'H', 'i', 2, 'j', 1),
    struct('name', 'GABAA_ss_rl', 'matrix', 'H', 'i', 8, 'j', 1),
    struct('name', 'GABAA_sp_sp', 'matrix', 'H', 'i', 2, 'j', 2),
    struct('name', 'GABAA_sp_si', 'matrix', 'H', 'i', 3, 'j', 2),
    struct('name', 'GABAA_sp_dp', 'matrix', 'H', 'i', 4, 'j', 2),
    struct('name', 'GABAA_si_si', 'matrix', 'H', 'i', 3, 'j', 3),
    struct('name', 'GABAA_si_sp', 'matrix', 'H', 'i', 2, 'j', 3),
    struct('name', 'GABAA_dp_dp', 'matrix', 'H', 'i', 4, 'j', 4),
    struct('name', 'GABAA_dp_di', 'matrix', 'H', 'i', 5, 'j', 4),
    struct('name', 'GABAA_dp_tp', 'matrix', 'H', 'i', 6, 'j', 4),
    struct('name', 'GABAA_di_di', 'matrix', 'H', 'i', 5, 'j', 5),
    struct('name', 'GABAA_di_dp', 'matrix', 'H', 'i', 4, 'j', 5),
    struct('name', 'GABAA_di_tp', 'matrix', 'H', 'i', 6, 'j', 5),
    struct('name', 'GABAA_tp_tp', 'matrix', 'H', 'i', 6, 'j', 6),
    struct('name', 'GABAA_tp_rl', 'matrix', 'H', 'i', 8, 'j', 6),
    struct('name', 'GABAA_rt_rt', 'matrix', 'H', 'i', 7, 'j', 7),
    struct('name', 'GABAA_rt_rl', 'matrix', 'H', 'i', 8, 'j', 7),
    struct('name', 'GABAA_rl_rl', 'matrix', 'H', 'i', 8, 'j', 8),
    struct('name', 'GABAA_rl_ss', 'matrix', 'H', 'i', 1, 'j', 8)
};

% GABA-B Parameters
GABAB_params = {
    struct('name', 'GABAB_ss_ss', 'matrix', 'Gb', 'i', 1, 'j', 1),
    struct('name', 'GABAB_ss_sp', 'matrix', 'Gb', 'i', 2, 'j', 1),
    struct('name', 'GABAB_ss_rl', 'matrix', 'Gb', 'i', 8, 'j', 1),
    struct('name', 'GABAB_sp_sp', 'matrix', 'Gb', 'i', 2, 'j', 2),
    struct('name', 'GABAB_si_si', 'matrix', 'Gb', 'i', 3, 'j', 3),
    struct('name', 'GABAB_si_sp', 'matrix', 'Gb', 'i', 2, 'j', 3),
    struct('name', 'GABAB_dp_dp', 'matrix', 'Gb', 'i', 4, 'j', 4),
    struct('name', 'GABAB_dp_di', 'matrix', 'Gb', 'i', 5, 'j', 4),
    struct('name', 'GABAB_dp_tp', 'matrix', 'Gb', 'i', 6, 'j', 4),
    struct('name', 'GABAB_tp_tp', 'matrix', 'Gb', 'i', 6, 'j', 6),
    struct('name', 'GABAB_tp_rl', 'matrix', 'Gb', 'i', 8, 'j', 6),
    struct('name', 'GABAB_rt_rt', 'matrix', 'Gb', 'i', 7, 'j', 7),
    struct('name', 'GABAB_rt_rl', 'matrix', 'Gb', 'i', 8, 'j', 7),
    struct('name', 'GABAB_rl_rl', 'matrix', 'Gb', 'i', 8, 'j', 8),
    struct('name', 'GABAB_rl_ss', 'matrix', 'Gb', 'i', 1, 'j', 8)
};

% NMDA Parameters
NMDA_params = {
    struct('name', 'NMDA_ss_sp', 'matrix', 'Hn', 'i', 2, 'j', 1),
    struct('name', 'NMDA_sp_sp', 'matrix', 'Hn', 'i', 2, 'j', 2),
    struct('name', 'NMDA_sp_si', 'matrix', 'Hn', 'i', 3, 'j', 2),
    struct('name', 'NMDA_sp_dp', 'matrix', 'Hn', 'i', 4, 'j', 2),
    struct('name', 'NMDA_si_sp', 'matrix', 'Hn', 'i', 2, 'j', 3),
    struct('name', 'NMDA_rl_ss', 'matrix', 'Hn', 'i', 1, 'j', 8)
};

neurotransmitter_names = {'AMPA', 'GABAA', 'GABAB', 'NMDA'};

%==========================================================================
% SECTION 3: SLEEP-STAGE SELECTION AND INPUT VALIDATION
%==========================================================================
sleep_stages = {'Wake'};  % Only process Wake stage
sleep_stage_patterns = {'*Wake*TCM*.mat'};  

% Display available files to verify path
fprintf('Looking for files in: %s\n', DCM_version);
all_files = dir(fullfile(DCM_version, '*.mat'));
fprintf('Found %d total .mat files\n', length(all_files));

% Check which stage files exist and include patient data.
fprintf('Checking which sleep stage files exist:\n');
valid_stages = {};
valid_indices = [];

for i = 1:length(sleep_stages)
    stage = sleep_stages{i};
    pattern = sleep_stage_patterns{i};
    files = dir(fullfile(DCM_version, pattern));
    
    % Count patient files only.
    patient_count = 0;
    for j = 1:length(files)
        filename = files(j).name;
        if is_patient_filename(filename)
            patient_count = patient_count + 1;
        end
    end
    
    fprintf('  %s: %d patient files matching "%s"\n', stage, patient_count, pattern);
    
    % A stage is valid when at least one patient DCM exists.
    if patient_count > 0
        valid_stages{end+1} = stage;
        valid_indices(end+1) = i;
    end
end

fprintf('Valid stages with patient data: %s\n', strjoin(valid_stages, ', '));

if isempty(valid_stages)
    error('No patient DCM files found for configured sleep-stage patterns.');
end

%==========================================================================
% SECTION 4: SIMULATION OPTIONS
%==========================================================================
mode = 'perc';   % 'abs' or 'perc'
vals = linspace(-1.0, 1.0, 5);
vals = [0, vals];

%==========================================================================
% SECTION 5: STAGE-WISE SENSITIVITY ANALYSIS
%==========================================================================
for idx = 1:length(valid_indices)
    stage_idx = valid_indices(idx);
    stage = sleep_stages{stage_idx};
    pattern = sleep_stage_patterns{stage_idx};
    fprintf('Processing Sleep Stage: %s\n', stage);
    
    %----------------------------------------------------------------------
    % 5.1 Load data for the current sleep stage
    %----------------------------------------------------------------------
    data_folder = fullfile(DCM_version);
    file_list   = dir(fullfile(data_folder, pattern));

    patient_DCM = {};  

    for file_idx = 1:length(file_list)
        filename  = file_list(file_idx).name;
        full_path = fullfile(data_folder, filename);
        is_patient_file = is_patient_filename(filename);
        if ~is_patient_file
            continue;
        end

        [loaded_ok, DCM] = load_dcm_and_psd(full_path, filename);
        if ~loaded_ok
            continue;
        end

        patient_DCM{end+1} = DCM;  
    end
    clear DCM

    % Verify data is loaded correctly
    n_patients = numel(patient_DCM);
    fprintf('  Loaded %d patient DCMs for stage %s\n', n_patients, stage);
    
    % Check if any patient data was found
    if isempty(patient_DCM)
        fprintf('  WARNING: No patient DCM data found for stage %s. Skipping this stage.\n', stage);
        continue;
    end

    if ~isfolder(atcm_dir)
        error('atcm_dir does not exist: %s', atcm_dir);
    end
    saved_working_dir = pwd;
    cleanup_atcm = onCleanup(@() cd(saved_working_dir));
    cd(atcm_dir);

    %----------------------------------------------------------------------
    % 5.2 Neurotransmitter-wise simulation sweeps
    %----------------------------------------------------------------------
    for group_idx = 1:numel(neurotransmitter_names)
        neurotransmitter_name = neurotransmitter_names{group_idx};
        switch neurotransmitter_name
            case 'AMPA'
                current_params = AMPA_params;
            case 'GABAA'
                current_params = GABAA_params;
            case 'GABAB'
                current_params = GABAB_params;
            case 'NMDA'
                current_params = NMDA_params;
            otherwise
                error('Unknown neurotransmitter group: %s', neurotransmitter_name);
        end

        fprintf('  Processing neurotransmitter group: %s (%d connections)\n', ...
            neurotransmitter_name, numel(current_params));

        for v = 1:numel(vals)
            current_val = vals(v);

            all_patient_sim_psd = [];
            patient_success_count = 0;
            patient_fail_count = 0;

            for patient_idx = 1:n_patients
                DCM = patient_DCM{patient_idx};
                current_psd = simulate_with_adjustment( ...
                    DCM, neurotransmitter_name, current_params, current_val, mode);

                if isempty(current_psd) || any(isnan(current_psd(:))) || any(isinf(current_psd(:)))
                    patient_fail_count = patient_fail_count + 1;
                    continue;
                end

                if isempty(all_patient_sim_psd)
                    all_patient_sim_psd = current_psd;
                elseif ~isequal(size(current_psd), size(all_patient_sim_psd))
                    patient_fail_count = patient_fail_count + 1;
                    continue;
                else
                    all_patient_sim_psd = all_patient_sim_psd + current_psd;
                end
                patient_success_count = patient_success_count + 1;
            end
            clear DCM

            if patient_fail_count > 0
                fprintf(['    %s adjustment %+.2f%%: %d/%d subjects OK' ...
                    ' (%d failed or size mismatch)\n'], ...
                    neurotransmitter_name, current_val * 100, ...
                    patient_success_count, n_patients, patient_fail_count);
            end

            if patient_success_count == 0 || isempty(all_patient_sim_psd)
                fprintf('    WARNING: No successful simulations for %s adjustment %.2f%% - skipping save.\n', ...
                    neurotransmitter_name, current_val * 100);
                continue;
            end

            sim_psd = all_patient_sim_psd / patient_success_count;

            if any(isnan(sim_psd(:))) || any(isinf(sim_psd(:)))
                fprintf('    WARNING: Invalid values in sim_psd for %s adjustment %.2f%%\n', ...
                    neurotransmitter_name, current_val * 100);
                continue;
            end

            psd_filename = fullfile(simulation_save_path, ...
                sprintf('%s_%s_Sim_Adjustment_%.2f.mat', ...
                stage, neurotransmitter_name, current_val * 100));
            parsave(psd_filename, sim_psd);
        end
    end
end

%==========================================================================
% SECTION 6: SAVE OUTPUTS
%==========================================================================
fprintf('Simulation sweeps completed.\n');