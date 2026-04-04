% PEB second-level: 22q11.2DS vs sibling. Runs spm_dcm_peb + BMC and saves GCM/PEB/BMA.

clearvars
close all

%==========================================================================
% SECTION 1: PATHS AND CONFIGURATION
%==========================================================================
scriptPath = fileparts(mfilename('fullpath'));
cd(scriptPath);
addpath(scriptPath);
addpath(fullfile(scriptPath, 'helperfunctions'));

base_path = '/Users/liobaberndt/Dropbox/UoE';
DCM_path  = fullfile(base_path, 'sleepdetectives', 'dcm', 'eLifeNick', 'DCM_v01');
PEB_path  = fullfile(base_path, 'sleepdetectives', 'dcm', 'eLifeNick', 'PEB_v08');
spm_path  = fullfile('/Users/liobaberndt/Desktop/spm_materials', 'spm12');

% { SPM field name, tag for filenames (stage_tag.mat) }
parameters = {
    {'Gb', 'GABAB'},
    {'H', 'GABAA'},
    {'Hn', 'NMDA'},
    {'H', 'AMPA'}
};

%==========================================================================
% SECTION 2: LOAD STAGE FILE LISTS
%==========================================================================
cd(DCM_path);
sleep_stages = {'Wake', 'N1', 'N2', 'N3', 'REM'};
stage_files = struct();
for stage = sleep_stages
    files = dir([stage{1} '*.mat']);
    stage_files.(stage{1}) = {files.name}';
end

%==========================================================================
% SECTION 3: PEB PER PARAMETER (per sleep stage × GABA/NMDA/AMPA field)
%==========================================================================
for stage_idx = 1:length(sleep_stages)
        stage_name = sleep_stages{stage_idx};
        if ~isfield(stage_files, stage_name)
            continue;
        end

    % Process each parameter type for this sleep stage
    for param_idx = 1:length(parameters)
        current_param = parameters{param_idx}{1};
        current_ver = parameters{param_idx}{2};

        Ver = [stage_name '_' current_ver];
        field = {current_param};

        % Files for this sleep stage (22q11.2DS vs sibling naming: pat / sib)
        all_files = stage_files.(stage_name);

        dcm_pat = {};
        dcm_sib = {};
        Ep_pat = [];
        Ep_sib = [];
        Cp_pat = [];
        Cp_sib = [];

        c = 1;  % 22q11.2DS (pat*)
        n = 1;  % siblings (sib*)

        for s = 1:length(all_files)
            filename = all_files{s};
            is_pat = ~isempty(regexp(filename, '^\w+-1(?:-\d+)?', 'once'));
            is_sib = ~is_pat;
            load(filename);

            if is_pat
                dcm_pat{c, 1} = filename;
                Ep_pat(c,:)  = full(spm_vec(DCM.Ep));
                Cp_pat(c,:)  = full(diag(DCM.Cp));
                c = c + 1;
            elseif is_sib
                dcm_sib{n, 1} = filename;
                Ep_sib(n,:)  = full(spm_vec(DCM.Ep));
                Cp_sib(n,:)  = full(diag(DCM.Cp));
                n = n + 1;
            end
        end

        Npat = c - 1;
        Nsib = n - 1;

        % Second level: sibling first (-1), then 22q11.2DS (+1)
        GCM = [dcm_sib; dcm_pat];
        M = struct();
        M.alpha = 1;
        M.beta  = 16;
        M.hE    = 0;
        M.hC    = 1/16;
        M.Q     = 'all';
        N = length(GCM);
        M.X = ones(N,2);
        M.Xnames = {'Mean', 'Group'};
        M.X(1:Nsib, 2) = -1;
        M.X((Nsib+1):(Nsib+Npat), 2) = 1;
        M.X(:,2:end) = M.X(:,2:end) - nanmean(M.X(:,2:end));
        M.X(:,2:end) = M.X(:,2:end)./nanstd(M.X(:,2:end),0);
        NaN_rows = any(isnan(M.X), 2);
        M.X(NaN_rows, :) = [];
        GCM(NaN_rows, :) = [];

        cd(PEB_path)
        save(['GCM_' Ver '.mat'],'GCM','M');
        cd(spm_path)
        [PEB,DCM] = spm_dcm_peb(GCM,M,field);
        cd(PEB_path)
        save(['PEB_' Ver '.mat'],'PEB');

        [BMC,PEB_best] = spm_dcm_bmc_peb(DCM,M,field);
        save(['PEBbest_' Ver '.mat'],'PEB_best');
        BMA = spm_dcm_peb_bmc(PEB_best);
        save(['BMA_' Ver '.mat'],'BMA')
        fprintf('  Done: %s\n', Ver);

end % param_idx

end % stage_idx