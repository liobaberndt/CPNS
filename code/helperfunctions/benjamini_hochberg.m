function q = benjamini_hochberg(p)
    q = NaN(size(p));
    valid = ~isnan(p);
    pv = p(valid);
    m = numel(pv);
    if m == 0
        return;
    end
    [ps, ord] = sort(pv(:));
    qs = ps .* m ./ (1:m)';
    for i = m-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end
    qs = min(qs, 1);
    qv = NaN(size(pv));
    qv(ord) = qs;
    q(valid) = qv;
end
