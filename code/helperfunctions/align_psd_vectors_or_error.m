function [a_out, b_out] = align_psd_vectors_or_error(a, b, context_file, label_a, label_b)
    if isempty(a) || isempty(b)
        error('Empty PSD (%s or %s) for %s.', label_a, label_b, context_file);
    end
    na = numel(a);
    nb = numel(b);
    if na ~= nb
        error('Length mismatch %s vs %s for %s (%d vs %d).', ...
            label_a, label_b, context_file, na, nb);
    end
    a_out = a(:).';
    b_out = b(:).';
end
