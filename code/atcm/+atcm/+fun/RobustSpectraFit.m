function [mnewspectra,Freqs,uncorrected,modelm,newspectra,FitPar ]=RobustSpectraFit(Freqs,singlespectra,PostSmooth)
% Robust average of trial spectra: fit log PSD vs log freq per trial (outliers down-weighted),
% then average residuals. Works on FULL frequency grid (no subset).
% Freqs: full frequency vector; singlespectra: trials x frequencies.

Nf   = length(Freqs);
N    = size(singlespectra,1);
% Avoid log(0): use frequencies > 0.1 Hz for the log-domain fit
f0   = find(Freqs >= 0.1, 1, 'first');
if isempty(f0), f0 = 1; end
f1   = Nf;
idx  = f0:f1;
Flog = log(max(Freqs(idx), 0.1));  % log-freq for fit

% Exclude 1.5–4.25 Hz from the linear fit (delta/theta bump) so slope is stable
fc1 = atcm.fun.findthenearest(1.5, Freqs(idx));
fc2 = atcm.fun.findthenearest(4.25, Freqs(idx));
if fc2 > fc1+1
    fitIdx = [1:fc1, fc2:length(idx)];
else
    fitIdx = 1:length(idx);
end
FitFreqs = Flog(fitIdx);

newspectra = zeros(N, Nf);
FitPar     = zeros(N, 2);
modelm     = zeros(N, Nf);

uncorrected = zeros(N, Nf);
for j = 1:N
    m = singlespectra(j,:);
    lm = log(max(m(idx), 1e-12));
    uncorrected(j, idx) = lm;
    if f0 > 1
        uncorrected(j, 1:f0-1) = log(max(m(1:f0-1), 1e-12));
    end
    if f1 < Nf
        uncorrected(j, f1+1:Nf) = log(max(m(f1+1:Nf), 1e-12));
    end
    Fitlm = lm(fitIdx);
    warning off;
    b = robustfit(FitFreqs, Fitlm);
    warning on;
    FitPar(j,:) = b;
    modelm_j = b(1) + b(2)*Flog;
    modelm(j, idx) = modelm_j;
    lm = lm - modelm_j;
    if PostSmooth > 0
        lm = atcm.fun.moving_average(lm, PostSmooth);
    end
    newspectra(j, idx) = lm;
    if f0 > 1
        newspectra(j, 1:f0-1) = uncorrected(j, 1:f0-1);
    end
    if f1 < Nf
        newspectra(j, f1+1:Nf) = uncorrected(j, f1+1:Nf);
    end
end
mnewspectra = mean(newspectra, 1);
Freqs = Freqs(:)';  % return original full grid
end

