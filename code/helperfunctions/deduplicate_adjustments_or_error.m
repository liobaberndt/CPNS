function entry = deduplicate_adjustments_or_error(entry, state, receptor)
    if isempty(entry.adjustments)
        return;
    end
    adj_key = round(entry.adjustments(:), 6);
    uniq = unique(adj_key);
    if numel(uniq) == numel(adj_key)
        return;
    end
    error('Duplicate adjustment values in %s/%s sweep (%d files, %d unique). Remove duplicate .mat files.', ...
        state, receptor, numel(adj_key), numel(uniq));
end
