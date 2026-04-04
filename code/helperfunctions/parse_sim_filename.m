function [ok, state, receptor, adjustment] = parse_sim_filename(file_name, states, receptor_types)
    ok = false;
    state = '';
    receptor = '';
    adjustment = NaN;

    [~, stem, ~] = fileparts(file_name);
    if contains(stem, '_Sim_Adjustment_') || contains(stem, '_Direct_Adjustment_')
        return;
    end
    expr = '^(Wake|N1|N2|N3|REM)_(AMPA|GABAA|GABAB|NMDA)_Adjustment_([-\d.]+)$';
    tokens = regexp(stem, expr, 'tokens', 'once');
    if isempty(tokens)
        return;
    end
    state_match = tokens{1};
    receptor_match = tokens{2};
    adjustment_token = tokens{3};

    state_ok = any(strcmp(states, state_match));
    if ~state_ok
        return;
    end

    receptor_ok = any(strcmp(receptor_types, receptor_match));
    if receptor_ok
        receptor = receptor_match;
    end
    if isempty(receptor)
        return;
    end

    adjustment = str2double(adjustment_token);
    if isnan(adjustment)
        return;
    end
    state = state_match;
    ok = true;
end
