function [Qp, Cp, Eh, F] = time_varying_optimization(DCM)
% Time-varying optimization using RK4 + Laplace at each time step
% This is a custom optimization that works with our time-varying transfer function

fprintf('Starting time-varying optimization...\n');

% Extract parameters
P = DCM.M.pE;
M = DCM.M;
U = DCM.xU;
Y = DCM.xY;

% Get frequencies from the data
freqs = M.Hz;
nf = length(freqs);

% Initialize optimization
max_iter = 50;
tol = 1e-6;
F_old = -inf;
F = -inf;

% Parameter bounds (simple bounds for now)
lb = -10 * ones(size(spm_vec(P)));
ub = 10 * ones(size(spm_vec(P)));

% Optimization loop
for iter = 1:max_iter
    fprintf('Iteration %d/%d\n', iter, max_iter);
    
    % Compute current transfer function
    try
        [Y_pred, ~, ~, ~, ~, ~] = atcm.fun.time_varying_tf(P, M, U);
        
        % Compute log-likelihood (simplified)
        if iscell(Y_pred) && ~isempty(Y_pred{1})
            Y_pred_vec = Y_pred{1};
            Y_obs_vec = Y.y{1};
            
            % Simple least squares objective
            residual = Y_obs_vec - Y_pred_vec;
            F = -sum(residual.^2);  % Negative because we want to maximize
            
            fprintf('  Current F = %.2e\n', F);
            
            % Check convergence
            if abs(F - F_old) < tol
                fprintf('Converged at iteration %d\n', iter);
                break;
            end
            F_old = F;
            
        else
            fprintf('  Warning: Invalid transfer function output\n');
            F = -1e10;  % Very bad objective
        end
        
    catch ME
        fprintf('  Error in transfer function: %s\n', ME.message);
        F = -1e10;  % Very bad objective
    end
    
    % Simple parameter update (gradient-free optimization)
    % This is a very basic approach - in practice you'd want proper gradients
    if iter < max_iter
        % Random walk in parameter space
        step_size = 0.1 * exp(-iter/20);  % Decreasing step size
        P_vec = spm_vec(P);
        P_vec = P_vec + step_size * randn(size(P_vec));
        P_vec = max(min(P_vec, ub), lb);  % Apply bounds
        P = spm_unvec(P_vec, P);
    end
end

% Return results
Qp = P;
Cp = eye(length(spm_vec(P)));  % Simple covariance
Eh = [];
F = F;

fprintf('Time-varying optimization complete. Final F = %.2e\n', F);

end

