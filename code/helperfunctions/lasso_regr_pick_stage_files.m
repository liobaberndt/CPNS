function files_out = lasso_regr_pick_stage_files(stage_all, subject_ids)
    files_out = cell(size(subject_ids));
    for i = 1:numel(subject_ids)
        sid = subject_ids{i};
        hit = '';
        for j = 1:numel(stage_all)
            p = dcm_wake_parse(stage_all{j});
            if p.ok && strcmp(p.token, sid)
                hit = stage_all{j};
                break;
            end
        end
        if isempty(hit)
            error('lasso_regr_pick_stage_files:MissingStageFile', ...
                'No file in stage list for subject id ''%s''.', sid);
        end
        files_out{i} = hit;
    end
end
