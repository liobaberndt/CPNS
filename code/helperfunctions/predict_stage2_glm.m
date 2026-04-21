function yhat = predict_stage2_glm(mdl, X)
    if isstruct(mdl) && isfield(mdl, 'model_class') && strcmp(mdl.model_class, 'custom_negbin_glm')
        if isempty(X)
            Z = ones(1, numel(mdl.beta));
        else
            Z = [ones(size(X, 1), 1), X];
        end
        yhat = exp(Z * mdl.beta);
        yhat = yhat(:);
    else
        yhat = predict(mdl, X);
        yhat = yhat(:);
    end
end
