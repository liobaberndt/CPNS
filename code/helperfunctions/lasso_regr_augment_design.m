function [X_pen, X_cov, names_pen, names_cov] = lasso_regr_augment_design(X, G, age_z, sex_v, param_names)
    [n, p] = size(X);
    if numel(G) ~= n || numel(age_z) ~= n || numel(sex_v) ~= n
        error('lasso_regr_augment_design:AugmentSize', ...
            'Covariates must match DCM rows (%d).', n);
    end
    if numel(param_names) ~= p
        error('lasso_regr_augment_design:AugmentNames', ...
            'param_names length (%d) ~= DCM columns (%d).', numel(param_names), p);
    end
    Xc = X - mean(X, 1, 'omitnan');
    Gcol = double(G(:));
    az = double(age_z(:));
    sv = double(sex_v(:));
    two_groups = numel(unique(Gcol)) > 1;

    names_cov = {};
    blk_cov = [];
    if two_groups
        blk_cov = [blk_cov, Gcol];
        names_cov{end+1} = 'group_22q';
    end
    blk_cov = [blk_cov, az, sv];
    names_cov{end+1} = 'age_z';
    names_cov{end+1} = 'sex';
    X_cov = blk_cov;
    names_cov = names_cov(:);

    names_pen = cell(p, 1);
    blk_pen = Xc;
    for j = 1:p
        names_pen{j} = param_names{j};
    end
    if two_groups
        Gi = Gcol .* Xc;
        blk_pen = [blk_pen, Gi];
        for j = 1:p
            names_pen{end+1} = sprintf('%s_x_group22q', param_names{j});
        end
    end
    X_pen = blk_pen;
    names_pen = names_pen(:);
end
