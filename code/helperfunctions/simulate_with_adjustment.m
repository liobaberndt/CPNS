function sim_psd = simulate_with_adjustment(DCM, neurotransmitter_name, param_list, adjustment, mode)
    sim_DCM = DCM;

    if strcmp(neurotransmitter_name, 'AMPA') || strcmp(neurotransmitter_name, 'GABAA')
        field_name = 'H';
    elseif strcmp(neurotransmitter_name, 'GABAB')
        field_name = 'Gb';
    elseif strcmp(neurotransmitter_name, 'NMDA')
        field_name = 'Hn';
    else
        fprintf('    WARNING: Unknown neurotransmitter type: %s\n', neurotransmitter_name);
        sim_psd = [];
        return;
    end

    if ~isfield(sim_DCM.Ep, field_name)
        fprintf('    WARNING: DCM.Ep does not contain field: %s\n', field_name);
        sim_psd = [];
        return;
    end

    M = sim_DCM.Ep.(field_name);
    if ~isnumeric(M) || ~ismatrix(M) || ~isequal(size(M), [8, 8])
        fprintf(['    WARNING: DCM.Ep.%s must be 8x8 numeric ' ...
            '(Wake pipeline); got %s of size %s.\n'], field_name, class(M), mat2str(size(M)));
        sim_psd = [];
        return;
    end
    sim_DCM.Ep.(field_name) = apply_nt_sweep(M, param_list, adjustment, mode);

    if ~isstruct(sim_DCM) || ~isfield(sim_DCM, 'Ep') || ~isstruct(sim_DCM.Ep) || ...
            ~isfield(sim_DCM, 'M') || ~isstruct(sim_DCM.M)
        sim_psd = [];
        return;
    end

    if ~isfield(sim_DCM.M, 'Hz') || isempty(sim_DCM.M.Hz)
        sim_DCM.M.Hz = 1:0.25:30;
    end

    if ~isfield(sim_DCM.M, 'nograph')
        sim_DCM.M.nograph = 1;
    end
    if ~isfield(sim_DCM.M, 'IS')
        sim_DCM.M.IS = @atcm.fun.alex_tf;
    elseif ischar(sim_DCM.M.IS) || isstring(sim_DCM.M.IS)
        sim_DCM.M.IS = str2func(char(sim_DCM.M.IS));
    end

    if ~isstruct(sim_DCM.xU)
        sim_DCM.xU = struct();
    end
    if ~isfield(sim_DCM.xU, 'X')
        sim_DCM.xU.X = [];
    end
    if ~isfield(sim_DCM.xU, 'dt')
        sim_DCM.xU.dt = 0.001;
    end

    try
        [Y, ~, ~, ~, ~, ~] = atcm.fun.alex_tf(sim_DCM.Ep, sim_DCM.M, sim_DCM.xU);
        sim_psd = Y;
        if iscell(sim_psd) && numel(sim_psd) == 1
            sim_psd = sim_psd{1};
        end

        if any(isnan(sim_psd(:))) || any(isinf(sim_psd(:))) || all(sim_psd(:) == 0)
            sim_psd = [];
        end
    catch ME
        fprintf('    ERROR in simulation: %s\n', ME.message);
        sim_psd = [];
    end
end
