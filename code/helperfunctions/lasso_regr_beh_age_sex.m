function [age_v, sex_v] = lasso_regr_beh_age_sex(beh)
    if ~isfield(beh, 'age')
        error('lasso_regr_beh_age_sex:MissingAge', 'behavioural data must include field ''age''.');
    end
    if ~isfield(beh, 'sex')
        error('lasso_regr_beh_age_sex:MissingSex', 'behavioural data must include field ''sex''.');
    end
    age_v = double(beh.age(:));
    raw = beh.sex;
    if isnumeric(raw)
        sex_v = double(raw(:));
        return;
    end
    sex_v = nan(numel(raw), 1);
    if iscell(raw)
        for ii = 1:numel(raw)
            sex_v(ii) = local_sex_char_to_num(raw{ii});
        end
    else
        for ii = 1:numel(raw)
            sex_v(ii) = local_sex_char_to_num(raw(ii));
        end
    end
end

function v = local_sex_char_to_num(s)
    if isnumeric(s)
        v = double(s(1));
        return;
    end
    t = lower(strtrim(char(string(s))));
    if isempty(t) || strcmp(t, 'nan')
        v = NaN;
    elseif any(strcmp(t, { 'f', 'female', 'woman', 'girl', '2' }))
        v = 1;
    elseif any(strcmp(t, { 'm', 'male', 'man', 'boy', '1' }))
        v = 0;
    else
        v = NaN;
    end
end
