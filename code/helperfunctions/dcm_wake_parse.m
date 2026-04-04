function p = dcm_wake_parse(filename)
    parts = regexp(filename, '_', 'split');
    if numel(parts) < 4
        p = struct('ok', false, 'token', '', 'is_patient', false);
        return;
    end
    tok = strtrim(lower(char(parts{4})));
    is_pat = ~isempty(regexp(tok, '^\w+-1(?:-\d+)?', 'once'));
    p = struct('ok', true, 'token', tok, 'is_patient', is_pat);
end
