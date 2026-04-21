function coefT = stage2_glm_coefficients(mdl)
    % Shared accessor keeps downstream code agnostic to model class.
    coefT = mdl.Coefficients;
end
