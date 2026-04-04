function beh_out = lasso_regr_beh_align_to_dcm(beh_in, W, n_sub)
    if ~isfield(beh_in, 'ids')
        error('lasso_regr_beh_align_to_dcm:MissingIds', ...
            'behavioural data must include field ''ids''.');
    end
    fields = fieldnames(beh_in);
    ids_norm = lasso_regr_ids_to_cellstr(beh_in.ids);
    dcm_ids = cell(n_sub, 1);
    for i = 1:n_sub
        parts = regexp(W{i}, '_', 'split');
        if numel(parts) < 4
            error('lasso_regr_beh_align_to_dcm:BadWakeName', ...
                'Wake filename needs subject id as 4th _-segment: %s', W{i});
        end
        dcm_ids{i} = strtrim(lower(char(parts{4})));
    end

    row_map = zeros(n_sub, 1);
    for i = 1:n_sub
        hit = find(strcmp(ids_norm, dcm_ids{i}), 1);
        if isempty(hit)
            error('lasso_regr_beh_align_to_dcm:BehIdMissing', ...
                'No behavioural row with ids matching ''%s'' (from %s).', dcm_ids{i}, W{i});
        end
        row_map(i) = hit;
    end
    if numel(unique(row_map)) ~= n_sub
        error('lasso_regr_beh_align_to_dcm:DuplicateIdMap', ...
            'Multiple DCM rows map to the same behavioural row; ids must be unique.');
    end

    beh_out = struct();
    for k = 1:numel(fields)
        f = fields{k};
        beh_out.(f) = beh_in.(f)(row_map);
    end
end

function c = lasso_regr_ids_to_cellstr(ids_raw)
    if iscell(ids_raw)
        c = cell(numel(ids_raw), 1);
        for ii = 1:numel(ids_raw)
            x = ids_raw{ii};
            if isnumeric(x)
                c{ii} = strtrim(lower(num2str(double(x(:)'))));
                c{ii} = strrep(c{ii}, '  ', ' ');
            else
                c{ii} = strtrim(lower(char(x)));
            end
        end
    elseif isnumeric(ids_raw)
        c = cell(numel(ids_raw), 1);
        for ii = 1:numel(ids_raw)
            c{ii} = strtrim(lower(num2str(double(ids_raw(ii)))));
        end
    elseif isa(ids_raw, 'string')
        c = cellstr(ids_raw(:));
        for ii = 1:numel(c)
            c{ii} = strtrim(lower(c{ii}));
        end
    else
        tmp = char(ids_raw);
        if size(tmp, 1) > 1
            c = cell(size(tmp, 1), 1);
            for ii = 1:size(tmp, 1)
                c{ii} = strtrim(lower(strtrim(tmp(ii, :))));
            end
        else
            c = {strtrim(lower(strtrim(tmp)))};
        end
    end
end
