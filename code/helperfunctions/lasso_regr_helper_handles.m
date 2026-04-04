function h = lasso_regr_helper_handles()
    names = { 'nested_cv_lasso_stage2', 'overdispersion_pearson_poisson', 'benjamini_hochberg' };
    h = struct();
    for k = 1:numel(names)
        nm = names{k};
        if isempty(which(nm))
            error('lasso_regr_helper_handles:MissingHelper', ...
                '%s not on path.', nm);
        end
        h.(nm) = str2func(nm);
    end
end
