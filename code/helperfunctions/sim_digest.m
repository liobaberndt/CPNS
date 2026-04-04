function d = sim_digest(v)
    k = min(10, numel(v));
    d = sum(v(1:k));
end
