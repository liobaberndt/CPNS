function coefT = stage2_glm_coefficients(mdl)
    if isstruct(mdl) && isfield(mdl, 'model_class') && strcmp(mdl.model_class, 'custom_negbin_glm')
        coefT = mdl.Coefficients;
    else
        coefT = mdl.Coefficients;
    end
end
