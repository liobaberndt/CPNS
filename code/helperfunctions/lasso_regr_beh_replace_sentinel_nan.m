function s = lasso_regr_beh_replace_sentinel_nan(s, sentinel)
    skipId = lower({ 'ids', 'id', 'subject_id', 'subjectid' });
    fn = fieldnames(s);
    for k = 1:numel(fn)
        f = fn{k};
        if any(strcmp(lower(f), skipId))
            continue;
        end
        v = s.(f);
        if isnumeric(v)
            v = double(v);
            bad = (v == sentinel) | (v == -double(sentinel));
            v(bad) = NaN;
            s.(f) = v;
        end
    end
end
