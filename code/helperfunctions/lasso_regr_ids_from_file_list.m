function ids = lasso_regr_ids_from_file_list(files)
    ids = cell(numel(files), 1);
    for i = 1:numel(files)
        p = dcm_wake_parse(files{i});
        if ~p.ok
            error('lasso_regr_ids_from_file_list:BadFile', ...
                'Expected subject id as 4th _-segment: %s', files{i});
        end
        ids{i} = p.token;
    end
end
