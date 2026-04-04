function code = lasso_regr_group_code(g)
    n = numel(g);
    code = zeros(n, 1);
    for i = 1:n
        if iscell(g)
            gi = char(string(g{i}));
        else
            gi = char(string(g(i)));
        end
        gi = lower(strtrim(gi));
        if isempty(gi)
            continue;
        end
        if contains(gi, '22q')
            code(i) = 1;
        elseif contains(gi, 'sib') || contains(gi, 'sibling') || contains(gi, 'control') ...
                || strcmp(gi, 'ctrl') || strcmp(gi, 'hc')
            code(i) = 2;
        end
    end
end
