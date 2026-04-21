% DCM parameters vs psychiatric/cognitive outcomes — nested CV LASSO, stage-2 GLM, BH-FDR.

function lasso_regr

clearvars

%==========================================================================
% SECTION 1: PATHS AND CONFIGURATION
%==========================================================================
scriptPath = fileparts(mfilename('fullpath'));
cd(scriptPath);
addpath(scriptPath);
addpath(fullfile(scriptPath, 'helperfunctions'));

base_path  = '/Users/liobaberndt/Dropbox/UoE';
model_path = fullfile(base_path, 'sleepdetectives', 'dcm', 'eLifeNick', 'DCM_v01');
fun_path   = fullfile(base_path, 'DCM_TCM', 'atcm');
psych_path = fullfile(base_path, 'sleepdetectives', 'data', 'behav_data');
addpath(fun_path);

fprintf('DCM path: %s\n', model_path);
fprintf('atcm path: %s\n', fun_path);
fprintf('Psych data path: %s\n', psych_path);

if isempty(which('nested_cv_lasso_stage2'))
    error('lasso_regr:MissingNestedCV', 'Add code/helperfunctions to the MATLAB path.');
end

%==========================================================================
% SECTION 2: DCM FILE LISTS (22q + SIBLINGS)
%==========================================================================
cd(model_path);

d = dir('Wake*.mat');
W_all = { d.name }';
d = dir('N1*.mat');
N1d_all = { d.name }';
d = dir('N2*.mat');
N2d_all = { d.name }';
d = dir('N3*.mat');
N3d_all = { d.name }';
d = dir('REM*.mat');
REMd_all = { d.name }';

W_pat = filter_patient_data(W_all);
W_sib = filter_sibling_data(W_all);
if isempty(W_pat)
    error('lasso_regr:NoPatientDcm', 'No patient Wake*.mat files under DCM path.');
end

ids_pat = lasso_regr_ids_from_file_list(W_pat);
ids_sib = lasso_regr_ids_from_file_list(W_sib);
ids_pat = sort(ids_pat);
ids_sib = sort(ids_sib);
if ~isempty(intersect(ids_pat, ids_sib))
    error('lasso_regr:DuplicateSubjectIds', 'Same id in patient and sibling Wake lists.');
end
subject_ids = [ids_pat(:); ids_sib(:)];
n_sub = numel(subject_ids);
n_pat = numel(ids_pat);
n_sib = numel(ids_sib);
fprintf('DCM subjects: n_22q=%d, n_sibling=%d, n_total=%d\n', n_pat, n_sib, n_sub);

W    = lasso_regr_pick_stage_files(W_all, subject_ids);
N1d  = lasso_regr_pick_stage_files(N1d_all, subject_ids);
N2d  = lasso_regr_pick_stage_files(N2d_all, subject_ids);
N3d  = lasso_regr_pick_stage_files(N3d_all, subject_ids);
REMd = lasso_regr_pick_stage_files(REMd_all, subject_ids);

%==========================================================================
% SECTION 3: LOAD DCM SPECTRA AND POSTERIOR PARAMETERS
%==========================================================================
for i = 1:n_sub
    load(fullfile(model_path, W{i}), 'DCM');
    Wake.model(i, :) = DCM.pred{:};
    Wake.data(i, :)  = DCM.xY.y{:};
    Wake.r2(i) = corr(squeeze(Wake.model(i, :)'), squeeze(Wake.data(i, :)')) .^ 2;
    Wake.params(i) = atcm.get_posteriors_tcm2024(DCM.Ep);

    load(fullfile(model_path, N1d{i}), 'DCM');
    N1.model(i, :) = DCM.pred{:};
    N1.data(i, :)  = DCM.xY.y{:};
    N1.r2(i) = corr(squeeze(N1.model(i, :)'), squeeze(N1.data(i, :)')) .^ 2;
    N1.params(i) = atcm.get_posteriors_tcm2024(DCM.Ep);

    load(fullfile(model_path, N2d{i}), 'DCM');
    N2.model(i, :) = DCM.pred{:};
    N2.data(i, :)  = DCM.xY.y{:};
    N2.r2(i) = corr(squeeze(N2.model(i, :)'), squeeze(N2.data(i, :)')) .^ 2;
    N2.params(i) = atcm.get_posteriors_tcm2024(DCM.Ep);

    load(fullfile(model_path, N3d{i}), 'DCM');
    N3.model(i, :) = DCM.pred{:};
    N3.data(i, :)  = DCM.xY.y{:};
    N3.r2(i) = corr(squeeze(N3.model(i, :)'), squeeze(N3.data(i, :)')) .^ 2;
    N3.params(i) = atcm.get_posteriors_tcm2024(DCM.Ep);

    load(fullfile(model_path, REMd{i}), 'DCM');
    REM.model(i, :) = DCM.pred{:};
    REM.data(i, :)  = DCM.xY.y{:};
    REM.r2(i) = corr(squeeze(REM.model(i, :)'), squeeze(REM.data(i, :)')) .^ 2;
    REM.params(i) = atcm.get_posteriors_tcm2024(DCM.Ep);
end

%==========================================================================
% SECTION 4: STACK PARAMETERS INTO MATRICES PER STAGE
%==========================================================================
[~, Wakemat, param_names_base] = param_struct_to_table(Wake.params);
[~, N1mat, ~] = param_struct_to_table(N1.params);
[~, N2mat, ~] = param_struct_to_table(N2.params);
[~, N3mat, ~] = param_struct_to_table(N3.params);
[~, REMmat, ~] = param_struct_to_table(REM.params);

%==========================================================================
% SECTION 5: BEHAVIOURAL DATA, ALIGN, FILTER, AUGMENT DESIGN
%==========================================================================
psych_mat = fullfile(psych_path, 'beh_psych_demo_data.mat');
if exist(psych_mat, 'file') ~= 2
    error('lasso_regr:MissingPsychMat', 'Expected %s', psych_mat);
end
psych_data = load(psych_mat);
if ~isfield(psych_data, 'data')
    error('lasso_regr:BadPsychMat', 'beh_psych_demo_data.mat must contain variable ''data''.');
end
beh_data_1 = psych_data.data;
req = { 'group', 'sleepprobsall', 'age', 'sex', 'ids' };
for ri = 1:numel(req)
    if ~isfield(beh_data_1, req{ri})
        error('lasso_regr:MissingBehField', 'beh_data must include field ''%s''.', req{ri});
    end
end

beh_data_1 = lasso_regr_beh_align_to_dcm(beh_data_1, W, n_sub);
beh_data_1 = lasso_regr_beh_replace_sentinel_nan(beh_data_1, 8888);

fields = fieldnames(beh_data_1);
sp = beh_data_1.sleepprobsall(:);
if ~isnumeric(sp)
    error('lasso_regr:SleepType', 'sleepprobsall must be numeric; got %s.', class(sp));
end
gc = lasso_regr_group_code(beh_data_1.group);
keep_mask = (gc > 0) & isfinite(sp);
keep_idx = find(keep_mask);
if isempty(keep_idx)
    error('lasso_regr:NoSubjects', 'No rows left after group + sleepprobsall filter (n_sub=%d).', n_sub);
end

[age_all, sex_all] = lasso_regr_beh_age_sex(beh_data_1);

beh_data = struct();
for ii = 1:numel(fields)
    field = fields{ii};
    v = beh_data_1.(field);
    if iscell(v)
        beh_data.(field) = v(keep_idx);
    elseif isnumeric(v)
        beh_data.(field) = v(keep_idx);
    else
        error('lasso_regr:BehFieldType', 'Unsupported type for field %s.', field);
    end
end

Wakemat = Wakemat(keep_idx, :);
N1mat   = N1mat(keep_idx, :);
N2mat   = N2mat(keep_idx, :);
N3mat   = N3mat(keep_idx, :);
REMmat  = REMmat(keep_idx, :);

group_22q = (gc == 1);
group_22q = group_22q(keep_idx);

age_kept = double(age_all(keep_idx));
sex_kept = double(sex_all(keep_idx));
age_z = (age_kept - mean(age_kept, 'omitnan')) ./ std(age_kept, 0, 'omitnan');
age_z(~isfinite(age_z)) = 0;

[Wakemat_pen, Wakemat_cov, names_pen, names_cov] = lasso_regr_augment_design(Wakemat, group_22q, age_z, sex_kept, param_names_base);
[N1mat_pen, N1mat_cov, ~, ~] = lasso_regr_augment_design(N1mat, group_22q, age_z, sex_kept, param_names_base);
[N2mat_pen, N2mat_cov, ~, ~] = lasso_regr_augment_design(N2mat, group_22q, age_z, sex_kept, param_names_base);
[N3mat_pen, N3mat_cov, ~, ~] = lasso_regr_augment_design(N3mat, group_22q, age_z, sex_kept, param_names_base);
[REMmat_pen, REMmat_cov, ~, ~] = lasso_regr_augment_design(REMmat, group_22q, age_z, sex_kept, param_names_base);

param_names = [names_cov; names_pen];
n_cov = numel(names_cov);

%==========================================================================
% SECTION 6: OUTCOME DEFINITIONS
%==========================================================================
outcomes = struct('field', {}, 'label', {}, 'lasso', {}, 'stage2', {}, 'y_xform', {});
o = 1;
outcomes(o) = struct('field', 'adhdsymsall', 'label', 'ADHD', 'lasso', 'poisson', 'stage2', 'negbin', 'y_xform', 'none'); o = o + 1;
outcomes(o) = struct('field', 'anyanxsymsall', 'label', 'Anxiety', 'lasso', 'poisson', 'stage2', 'negbin', 'y_xform', 'none'); o = o + 1;
outcomes(o) = struct('field', 'asqtotalsymsall', 'label', 'ASD', 'lasso', 'poisson', 'stage2', 'negbin', 'y_xform', 'none'); o = o + 1;
outcomes(o) = struct('field', 'sleepprobsall', 'label', 'Sleep_Problems', 'lasso', 'poisson', 'stage2', 'negbin', 'y_xform', 'none'); o = o + 1;
outcomes(o) = struct('field', 'total_pe', 'label', 'Psychotic_Experiences', 'lasso', 'binomial', 'stage2', 'logistic', 'y_xform', 'none'); o = o + 1;
outcomes(o) = struct('field', 'fsiqall', 'label', 'Full_Spectrum_IQ', 'lasso', 'gaussian', 'stage2', 'linear', 'y_xform', 'none'); o = o + 1;

extra_outcomes = {
    'PerseverativeErrorsStandardall',  'Perseverative_Errors',  'poisson',  'negbin',   'none'
    'NonperseverativeErrorsStall',      'Nonperseverative_Err',  'poisson',  'negbin',   'none'
    'RVPAstandardall',                  'RVP_Aprime',            'gaussian', 'linear',   'none'
    'SOCprobsSTANDall',                 'SOC_Problems_Solved',   'gaussian', 'linear',   'none'
    'MTStotalall',                      'Movement_Time_MTS',     'gaussian', 'linear',   'log1p'
    'RTIreactSTANall',                  'RTI_Reaction_ProcSpd',  'gaussian', 'linear',   'log1p'
    'RTImoveSTANall',                   'RTI_Movement',          'gaussian', 'linear',   'log1p'
    'SOCinitialSTANDall',               'SOC_Initial_Time',      'gaussian', 'linear',   'log1p'
    'SOCsubSTANDall',                   'SOC_Subsequent_Time',   'gaussian', 'linear',   'log1p'
    };
for r = 1:size(extra_outcomes, 1)
    if isfield(beh_data, extra_outcomes{r, 1})
        outcomes(end+1) = struct('field', extra_outcomes{r, 1}, 'label', extra_outcomes{r, 2}, ...
            'lasso', extra_outcomes{r, 3}, 'stage2', extra_outcomes{r, 4}, 'y_xform', extra_outcomes{r, 5});
    end
end

stage_mats_pen = { Wakemat_pen, N1mat_pen, N2mat_pen, N3mat_pen, REMmat_pen };
stage_mats_cov = { Wakemat_cov, N1mat_cov, N2mat_cov, N3mat_cov, REMmat_cov };
stage_names = { 'Wake', 'N1', 'N2', 'N3', 'REM' };

psych_outcome_labels = { 'ADHD', 'Anxiety', 'ASD', 'Sleep_Problems', 'Psychotic_Experiences' };

%==========================================================================
% SECTION 7: NESTED CV LASSO, STABILITY SELECTION, STAGE-2 GLM
%==========================================================================
helper_f = lasso_regr_helper_handles();

regression_results = struct();
coef_table_rows = {};

outerK = 5;
innerK = 5;

for oi = 1:numel(outcomes)
    oc = outcomes(oi);
    if ~isfield(beh_data, oc.field)
        continue;
    end

    y_raw = beh_data.(oc.field);
    if isnumeric(y_raw)
        y = double(y_raw(:));
    else
        y = double(cell2mat(y_raw(:)));
        y = y(:);
    end

    nDesign = size(stage_mats_pen{1}, 1);
    if numel(y) ~= nDesign
        continue;
    end

    if strcmp(oc.stage2, 'logistic')
        if any(y ~= 0 & y ~= 1)
            y = double(y > 0);
        end
    elseif strcmp(oc.lasso, 'poisson') || strcmp(oc.stage2, 'negbin')
        y = max(round(y), 0);
    end

    if strcmp(oc.y_xform, 'log1p')
        y_model = real(log1p(max(y, -1 + eps)));
        y_model(~isfinite(y_model)) = NaN;
    else
        y_model = y;
    end

    if any(cellfun(@(s) strcmp(oc.field, s), extra_outcomes(:, 1)))
        y_model(abs(y_model) >= 8000) = NaN;
    end

    fprintf('-- Outcome %s (%s; LASSO=%s, stage2=%s) --\n', oc.label, oc.field, oc.lasso, oc.stage2);

    for si = 1:numel(stage_names)
        Xpen = stage_mats_pen{si};
        Xcov = stage_mats_cov{si};
        try
            res = helper_f.nested_cv_lasso_stage2(Xpen, y_model, oc.lasso, oc.stage2, ...
                'outerK', outerK, 'innerK', innerK, 'XCov', Xcov);
        catch ME
            warning('lasso_regr:NestedCVFailed', '%s / %s: %s', oc.label, stage_names{si}, ME.message);
            continue;
        end

        regression_results.(oc.label).(stage_names{si}) = res;
        fprintf('  %s: nested CV r=%.3f, mean|err|=%.4g, n_stable_pen=%d (>= %d/%d folds), n_cov=%d\n', ...
            stage_names{si}, res.cv_r, res.cv_mae_mean, sum(res.active_pen), res.min_votes, outerK, res.n_cov);

        mdl = res.final_mdl;
        if ~isempty(mdl)
            coefT = stage2_glm_coefficients(mdl);
            pidx = res.pred_idx_final(:);
            for rr = 1:numel(pidx)
                if rr + 1 > height(coefT)
                    break;
                end
                coef_table_rows(end+1, :) = { oc.label, stage_names{si}, param_names{pidx(rr)}, coefT.Estimate(rr + 1), coefT.pValue(rr + 1) };
            end
        end
    end
end

%==========================================================================
% SECTION 8: BH-FDR AND FIGURES
%==========================================================================
adj_coef = coef_table_rows;
if ~isempty(adj_coef)
    for ri = 1:size(adj_coef, 1)
        adj_coef{ri, 6} = NaN;
    end
end

if ~isempty(adj_coef)
    psych_mask = cellfun(@(lab) any(strcmp(lab, psych_outcome_labels)), adj_coef(:, 1));
    for si = 1:numel(stage_names)
        stg = stage_names{si};
        for domain = 1:2
            if domain == 1
                mask = strcmp(adj_coef(:, 2), stg) & psych_mask;
            else
                mask = strcmp(adj_coef(:, 2), stg) & ~psych_mask;
            end
            idx = find(mask);
            if isempty(idx)
                continue;
            end
            pv = cell2mat(adj_coef(idx, 5));
            qv = helper_f.benjamini_hochberg(pv);
            for k = 1:numel(idx)
                adj_coef{idx(k), 6} = qv(k);
            end
        end
    end
end

for ri = 1:size(adj_coef, 1)
    if size(adj_coef, 2) < 6 || isnan(adj_coef{ri, 6}) || adj_coef{ri, 6} >= 0.05
        continue;
    end
    lab = adj_coef{ri, 1};
    stg = adj_coef{ri, 2};
    pname = adj_coef{ri, 3};
    j = find(strcmp(param_names, pname), 1);
    if isempty(j)
        continue;
    end
    oc_list = outcomes(strcmp({ outcomes.label }, lab));
    if numel(oc_list) ~= 1
        continue;
    end
    oc = oc_list(1);
    yplt = beh_data.(oc.field);
    if isnumeric(yplt)
        yplt = double(yplt(:));
    else
        yplt = double(cell2mat(yplt(:)));
    end
    if strcmp(oc.stage2, 'logistic') && any(yplt ~= 0 & yplt ~= 1)
        yplt = double(yplt > 0);
    end
    if strcmp(oc.y_xform, 'log1p')
        yplt = real(log1p(max(double(yplt(:)), -1 + eps)));
        yplt(~isfinite(yplt)) = NaN;
    end
    xv = [];
    for si = 1:numel(stage_names)
        if strcmp(stage_names{si}, stg)
            if j <= n_cov
                xv = stage_mats_cov{si}(:, j);
            else
                xv = stage_mats_pen{si}(:, j - n_cov);
            end
            break;
        end
    end
    if isempty(xv)
        continue;
    end
    figure;
    scatter(xv, yplt, 36, 'filled');
    xlabel(param_names{j});
    ylabel(lab);
    title(sprintf('%s | %s vs %s (q=%.3g)', stg, param_names{j}, lab, adj_coef{ri, 6}));
    grid on;
end

end
