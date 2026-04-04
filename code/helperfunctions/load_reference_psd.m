function psd = load_reference_psd(file_path, candidate_fields)
    loaded = load(file_path);
    for i = 1:length(candidate_fields)
        field = candidate_fields{i};
        if isfield(loaded, field)
            psd = normalize_psd_vector(loaded.(field), sprintf('%s (%s)', file_path, field));
            return;
        end
    end
    error('In %s, expected one of: %s', file_path, strjoin(candidate_fields, ', '));
end
