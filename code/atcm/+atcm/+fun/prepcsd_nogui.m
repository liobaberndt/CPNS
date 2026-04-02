function DCM = prepcsd_nogui(DCM)
% Non-GUI version of prepcsd for cluster computing
% This is a modified version of atcm.fun.prepcsd that removes GUI dependencies
% for cluster computing environments

% Set defaults and Get D filename
%-------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    error('Please specify data and trials');
end

% ensure spatial modes have been computed
%-------------------------------------------------------------------------
try
    DCM.M.U;
catch
    error('Please estimate this model first');
end

% Add SPM path directly
spm_path = '/lustre/projects/Research_Project-T122332/spm_materials';
if ~exist('spm_eeg_load', 'file')
    fprintf('Adding SPM from: %s\n', spm_path);
    if exist(spm_path, 'dir')
        addpath(genpath(spm_path));
        try
            spm('defaults', 'eeg');
            fprintf('Successfully initialized SPM with EEG defaults\n');
        catch spm_err
            fprintf('Found SPM directory but failed to initialize: %s\n', spm_err.message);
        end
    else
        fprintf('WARNING: SPM directory not found at: %s\n', spm_path);
    end
end

% load D - without GUI fallbacks
%--------------------------------------------------------------------------
try
    D = spm_eeg_load(Dfile);
catch loadErr1
    try
        [p,f] = fileparts(Dfile);
        D = spm_eeg_load(f);
        DCM.xY.Dfile = fullfile(pwd,f);
    catch loadErr2
        % Try loading from the specified project path
        try
            projectPath = '/lustre/projects/Research_Project-T122332/UoE/sleepdetectives/data/elifeNick';
            [~,f,e] = fileparts(Dfile);
            fullPath = fullfile(projectPath, [f e]);
            D = spm_eeg_load(fullPath);
            DCM.xY.Dfile = fullPath;
        catch loadErr3
            % Try one more approach - check if we need to add .mat extension
            try
                if ~endsWith(Dfile, '.mat')
                    fullPath = [Dfile '.mat'];
                else
                    [~,f,~] = fileparts(Dfile);
                    fullPath = fullfile(projectPath, [f '.mat']);
                end
                D = spm_eeg_load(fullPath);
                DCM.xY.Dfile = fullPath;
            catch loadErr4
                error('Could not load dataset: %s. Please check file path.', Dfile);
            end
        end
    end
end

% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
if ~isfield(DCM.xY, 'modality')
    [mod, list] = modality(D, 0, 1);

    if isequal(mod, 'Multimodal')
        % Non-GUI alternative: just use the first modality
        DCM.xY.modality = list{1};
        fprintf('Multiple modalities found. Using: %s\n', DCM.xY.modality);
    else
        DCM.xY.modality = mod;
    end
end

if ~isfield(DCM.xY, 'Ic')
    DCM.xY.Ic = setdiff(D.meegchannels(DCM.xY.modality), D.badchannels);
else
    fprintf('Channel subspace selection: Using channel %d\n', DCM.xY.Ic);
end

Ic = DCM.xY.Ic;
Nc = length(Ic);
Nm = size(DCM.M.U,2);
DCM.xY.Ic = Ic;

% options
%--------------------------------------------------------------------------
try
    DT = DCM.options.D;
catch
    DT = 1;
end
try
    trial = DCM.options.trials;
catch
    trial = D.nconditions;
end

% alex - allow [] to specify all
if isempty(DCM.options.trials)
    trial = 1;
end

% get peristimulus times
%--------------------------------------------------------------------------
DCM.xY.Time = 1000*D.time; % ms
T1 = DCM.options.Tdcm(1);
T2 = DCM.options.Tdcm(2);

T1 = atcm.fun.findthenearest(T1, DCM.xY.Time);
T2 = atcm.fun.findthenearest(T2, DCM.xY.Time);

% Check for baseline specification
DOBASE = 0;
try
    baseT1 = DCM.options.baseTdcm(1);
    baseT2 = DCM.options.baseTdcm(2);
    [i, baseT1] = min(abs(DCM.xY.Time - baseT1));
    [i, baseT2] = min(abs(DCM.xY.Time - baseT2));
    baseIt = [baseT1:DT:baseT2];
    DCM.xY.baseIt = baseIt;
    DOBASE = 1;
end

% Time [ms] of down-sampled data
%----------------------------------------------------------------------
It = [T1:DT:T2]';               % indices - bins
DCM.xY.pst = DCM.xY.Time(It);   % PST
DCM.xY.It = It;                 % Indices of time bins
DCM.xY.dt = DT/D.fsample;       % sampling in seconds
Nb = length(It);                % number of bins
fprintf('Time Indices=%d:%d in steps of %d\n', T1, T2, DT);

% get frequency range
%--------------------------------------------------------------------------
Hz1 = DCM.options.Fdcm(1);      % lower frequency
Hz2 = DCM.options.Fdcm(2);      % upper frequency

% Krish Hack to allow finer frequency resolution
Fstep = 1;
try
    Fstep = DCM.options.FrequencyStep;
    fprintf('Using fstep=%f\n', Fstep);
catch
    Fstep = 1;
end

% Frequencies
%--------------------------------------------------------------------------
DCM.xY.Hz = (Hz1:Fstep:Hz2);    % Frequencies
Nf = length(DCM.xY.Hz);         % number of frequencies
Ne = length(trial);             % number of ERPs

% Cross spectral density for each trial type
%==========================================================================
condlabels = D.condlist;

try
    DCM.xY.y = DCM.Kinitial.ForceSpectra;
    fprintf('Using user-specified spectra...\n');
    DCM.xY.U = DCM.M.U;
    DCM.xY.code = condlabels(trial);
    return;
end

try 
    bpf = DCM.options.UseButterband;
end

% Initialize storage for actual conditions (for window approach)
actual_conditions = cell(1, Ne);

for i = 1:Ne
    clear Pfull
    
    cond_name = condlabels{trial(i)};
    if iscell(cond_name)
        cond_name = cond_name{1};
    end
    
    % trial indices
    %----------------------------------------------------------------------
    if isfield(DCM.options, 'window') && DCM.options.window
        % Window approach: TrialIndices are direct absolute trial indices
        if isfield(DCM.options, 'TrialIndices')
            c = DCM.options.TrialIndices;  % Direct trial indices
            % Validate trial indices are within valid range
            valid_trials = 1:size(D,3);
            invalid = c < 1 | c > size(D,3);
            if any(invalid)
                error('Invalid trial indices: %s. Valid range is 1:%d', mat2str(c(invalid)), size(D,3));
            end
            
            % Determine actual condition from selected trials and update metadata
            % If knownStage is provided (from grouping logic), use it instead of reading from D
            actual_cond = '';
            if isfield(DCM.options, 'knownStage') && ~isempty(DCM.options.knownStage)
                actual_cond = DCM.options.knownStage;
                cond_name = actual_cond;
                actual_conditions{i} = actual_cond;
            else
                % Fallback: read from D.trials if knownStage not provided
                try
                    if isfield(D, 'trials') && ~isempty(D.trials)
                        trial_conditions = cell(length(c), 1);
                        for idx = 1:length(c)
                            if c(idx) <= length(D.trials) && c(idx) > 0
                                trial_conditions{idx} = D.trials(c(idx)).label;
                            end
                        end
                        unique_conds = unique(trial_conditions);
                        if length(unique_conds) == 1
                            actual_cond = unique_conds{1};
                            cond_name = actual_cond;
                            actual_conditions{i} = actual_cond;
                        else
                            actual_conditions{i} = unique_conds{1};  % Use first condition if mixed
                        end
                    end
                catch
                    % Skip if trial condition info not available
                    actual_conditions{i} = cond_name;
                end
            end
        else
            error('Window approach requires TrialIndices to be set');
        end
    else
        actual_conditions{i} = cond_name;  % Standard approach uses original condition
        % Standard approach: filter by condition first, then index
        if ~isempty(DCM.options.trials)
            c = D.indtrial(condlabels(trial(i)), 'GOOD');
        else
            c = 1:size(D,3); 
        end
        
        if isfield(DCM.options, 'TrialIndices')
            c = c(DCM.options.TrialIndices);
        end
    end
    
    Nt = length(c);
    
    % Get data
    %----------------------------------------------------------------------
    P = zeros(Nf, Nm, Nm);
    if DOBASE == 1
        Pbase = zeros(Nf, Nm, Nm);
    end
    Fs = 1000/(DCM.xY.Time(2) - DCM.xY.Time(1));
    
    for j = 1:Nt
        % Get Y for this condition
        Y = full(double(D(Ic, It, c(j))'));
        Ymod = Y;
        
        if DOBASE == 1
            Ybase = full(double(D(Ic, baseIt, c(j))' * DCM.M.U));
        end
        
        try
            DoButterband = DCM.options.UseButterband;
            
            for ch = 1:size(Y, 2)
                Y(:, ch) = atcm.fun.bandpassfilter(Y(:, ch), 1/DCM.xY.dt, DoButterband)';
            end 
            
            if DOBASE
                for ch = 1:size(Ybase, 2)
                    Ybase(:, ch) = atcm.fun.bandpassfilter(Ybase(:, ch), 1/DCM.xY.dt, DoButterband)';
                end
            end
        end
        
        Ymod = Y;
        series{i}(j, :, :) = Y;
        
        UseWelch = 0;
        try 
            UseWelch = DCM.options.UseWelch;
        catch
            UseWelch = 0;
        end
        
        if (UseWelch == 1010)
            if isfield(DCM.options, 'vmd') && DCM.options.vmd
                [Ymod] = vmd(real(Ymod), 'NumIMFs', 5);
                
                [uu, ~, ~] = spm_svd(cov(Ymod'));                                        
                pcc = uu' * Ymod;
                nn = min(12, size(pcc, 1));
                pcc = pcc(1:nn, :);
                
                % autoregressive spectral method for VMD components
                for ii = 1:nn
                    Pfc(ii, :) = pyulear(pcc(ii, :), 2, DCM.xY.Hz, 1/DCM.xY.dt);
                end
                
                Pf(:, 1, 1) = sum(Pfc, 1);
                Pfull(j, :, :, :) = full(Pf);
            else
                % Simple fft(x)
                FFTSmooth = 0;
                try
                    FFTSmooth = DCM.options.FFTSmooth;
                end
                
                if UseWelch == 1010
                    % Alex fft
                    SpecFun = DCM.options.SpecFun;
                    
                    if strcmp(char(SpecFun), 'atcm.fun.AfftSmooth') && isfield(DCM.options, 'FFTbins')
                        [Pf, F] = SpecFun(Ymod', 1/DCM.xY.dt, DCM.xY.Hz, DCM.options.FFTbins);
                    else
                        [Pf, F] = SpecFun(Ymod', 1/DCM.xY.dt, DCM.xY.Hz);
                        
                        if j == 1
                            FMAT = atcm.fun.asinespectrum(DCM.xY.Hz, D.time(It) * 1000);
                        end
                    end
                    
                    if size(Pf, 1) == 1
                        Pf = Pf';
                    end
                    
                    Pfull(j, :, :, :) = full(Pf);
                end
            end
            
            % the base period spectrum
            if DOBASE
                [Pfb, F] = atcm.fun.Afft(Ybase', 1/DCM.xY.dt, DCM.xY.Hz);
                if FFTSmooth > 0
                    [dPfb, in] = atcm.fun.padtimeseries(Pfb);
                    dPfb = atcm.fun.HighResMeanFilt(dPfb, 1, FFTSmooth);
                    Pfb = dPfb(in);
                end
                Pfullbase(j, :, :, :) = full(Pfb);
            end
        end
    end
    
    % average trials for this condition
    if (UseWelch == 1010)
        if isfield(DCM.options, 'BeRobust') && DCM.options.BeRobust
            fprintf('Robust fitting\n');
            [mnewspectra, fq, unc, ~, ~, FitPar] = atcm.fun.RobustSpectraFit(F, Pfull, 2);
            mnewspectra = exp(mnewspectra);
            mnewspectra = mnewspectra - min(mnewspectra);
            Pfull = mnewspectra;
            
            DCM.xY.trial_spectra = unc;
            
            if DOBASE
                Pfullbase = squeeze(spm_robust_average(Pfullbase));
            end
        else
            DCM.xY.trial_spectra = Pfull;
            
            if size(Pfull, 1) > 1 && Nm == 1
                % retain first eigenmode
                m = 1;
                [u, s, v] = spm_svd(Pfull', 1);
                Pfull = u(:, m) * s(m, m) * mean(v(:, m));
                Pfull = full(Pfull)';
            elseif size(Pfull, 1) > 1 && Nm > 1
                % First try robust averaging with warning capture
                warning_state = warning('query', 'all');
                warning('off', 'all'); % Turn off warnings temporarily
                try
                    Pfull_robust = spm_robust_average(Pfull);
                    Pfull = squeeze(Pfull_robust);
                catch
                    % Fall back to regular averaging if robust fails
                    Pfull = squeeze(mean(Pfull, 1));
                end
                warning(warning_state); % Restore warning state
            else
                Pfull = squeeze(Pfull);
            end
            
            if DOBASE
                Pfullbase = squeeze(spm_robust_average(Pfullbase));
            end
        end
        
        if FFTSmooth > 0
            if Nm == 1
                Pfull = Pfull';
            end
            for nchanx = 1:size(Pf, 2)
                for nchany = 1:size(Pf, 3)
                    Pfxy = Pfull(:, nchanx, nchany);
                    [dPf, in] = atcm.fun.padtimeseries(Pfxy);
                    dPf = atcm.fun.HighResMeanFilt(dPf, 1, FFTSmooth);
                    Pfxy = dPf(in);
                    Pfull(:, nchanx, nchany) = Pfxy;
                end
            end
            if Nm == 1
                Pfull = Pfull';
            end
        end
        
        if DOBASE == 1
            Pfull = Pfull - Pfullbase;
        end
        
        P = Pfull;
        
        P(isinf(P)) = 0;
        P(isnan(P)) = 0;
        P = P - min(P(:));
        
        if isfield(DCM.options, 'han') && DCM.options.han
            % apply hanning window
            if Nc == 1
                % Use built-in rescale if available, otherwise use a simple version
                if exist('rescale', 'file')
                    P = P .* rescale(kaiser(Nf, 2.5), .01, 1)'.^.2;
                else
                    wind = kaiser(Nf, 2.5);
                    wind = 0.01 + (1 - 0.01) * (wind - min(wind)) / (max(wind) - min(wind));
                    P = P .* wind'.^.2;
                end
            else
                if exist('rescale', 'file')
                    wind = rescale(kaiser(Nf, 2.5), .01, 1)';
                else
                    wind = kaiser(Nf, 2.5);
                    wind = 0.01 + (1 - 0.01) * (wind - min(wind)) / (max(wind) - min(wind));
                    wind = wind';
                end
                P = P .* repmat(wind(:), [1 Nc, Nc]).^.2;
            end
        end
        
        try
            DCM.xY.csd{i} = P';
        catch
            DCM.xY.csd{i} = P;
        end
    end
end

% place cross-spectral density in xY.y
%==========================================================================
try
    DCM.xY.y = spm_cond_units(DCM.xY.csd, 'csd');
catch
    % Fallback if spm_cond_units is not available
    DCM.xY.y = DCM.xY.csd;
end

% KRISH HACK! If No normalisation is selected, do not normalise by power
DoNormalise = 1;
try
    DoNormalise = DCM.options.DoNormalise;
catch
    DoNormalise = 1;
end

if DoNormalise == 0
    try
        DCM.xY.y = krish_cond_unitsNONORMALISE(DCM.xY.csd, 'csd');
    catch
        DCM.xY.y = DCM.xY.csd;
    end
end

if DoNormalise == 2
    DCM.xY.y = DCM.xY.csd;
end

% output the series
DCM.xY.series = series;

DCM.xY.U = DCM.M.U;
% Use actual conditions if window approach was used, otherwise use original
if isfield(DCM.options, 'window') && DCM.options.window && exist('actual_conditions', 'var')
    DCM.xY.code = actual_conditions;
else
    DCM.xY.code = condlabels(trial);
end

try
    DCM.xY.FMAT = FMAT;
end

end 