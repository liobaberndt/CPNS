function psd = normalize_psd_vector(psd_raw, context_label)
    if nargin < 2 || isempty(context_label)
        context_label = 'sim_psd';
    end
    if iscell(psd_raw)
        if numel(psd_raw) ~= 1
            error('%s: sim_psd as cell must have exactly one element.', context_label);
        end
        psd_raw = psd_raw{1};
    end
    if isempty(psd_raw) || ~isnumeric(psd_raw)
        error('%s: expected non-empty numeric vector.', context_label);
    end
    if ~isreal(psd_raw)
        error('%s: spectrum must be real-valued.', context_label);
    end
    if ndims(psd_raw) > 2
        error('%s: spectrum must be 1D (ndims=%d).', context_label, ndims(psd_raw));
    end
    if size(psd_raw, 1) > 1 && size(psd_raw, 2) > 1
        error('%s: spectrum must be a vector, not a %d×%d matrix.', ...
            context_label, size(psd_raw, 1), size(psd_raw, 2));
    end
    psd = psd_raw(:).';
end
