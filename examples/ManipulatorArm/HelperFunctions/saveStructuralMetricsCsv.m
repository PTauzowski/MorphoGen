function saveStructuralMetricsCsv(metrics_ref, metrics_final, configs, volFrac, resultRoot)
% SAVESTRUCTURALMETRICSCSVS  Write structural performance metrics to CSV.
%
%   saveStructuralMetricsCsv(metrics_ref, metrics_final, configs, volFrac, resultRoot)
%
%   Saves one row per load configuration plus an aggregate MAX row.
%   Ratios are final/reference.
%
%   Inputs:
%     metrics_ref    struct from evaluateStructuralPerformance (x = ones, penal=1)
%     metrics_final  struct from evaluateStructuralPerformance (x = x_final)
%     configs        {nConfigs x 1} cell of config structs with .name and .label
%     volFrac        scalar  final volume fraction (mean of x_final)
%     resultRoot     char/string  output directory

    nConfigs = numel(configs);
    rows     = cell(nConfigs + 1, 1);

    for k = 1:nConfigs
        row.config_name  = string(configs{k}.name);
        row.vol_frac     = volFrac;
        row.sHM_ref      = metrics_ref.sHM_perConfig(k);
        row.sHM_final    = metrics_final.sHM_perConfig(k);
        row.sHM_ratio    = metrics_final.sHM_perConfig(k) / max(metrics_ref.sHM_perConfig(k), eps);
        row.u_ref        = metrics_ref.u_perConfig(k);
        row.u_final      = metrics_final.u_perConfig(k);
        row.u_ratio      = metrics_final.u_perConfig(k) / max(metrics_ref.u_perConfig(k), eps);
        rows{k} = row;
    end

    % Aggregate MAX row
    row.config_name  = "ALL_MAX";
    row.vol_frac     = volFrac;
    row.sHM_ref      = metrics_ref.sHM_max;
    row.sHM_final    = metrics_final.sHM_max;
    row.sHM_ratio    = metrics_final.sHM_max / max(metrics_ref.sHM_max, eps);
    row.u_ref        = metrics_ref.u_max;
    row.u_final      = metrics_final.u_max;
    row.u_ratio      = metrics_final.u_max / max(metrics_ref.u_max, eps);
    rows{end} = row;

    T = struct2table(vertcat(rows{:}));
    writetable(T, fullfile(resultRoot, 'structural_metrics.csv'));
    fprintf('[structural metrics] sHM: %.3e -> %.3e (x%.2f)  u: %.3e -> %.3e (x%.2f)  vf=%.3f\n', ...
        metrics_ref.sHM_max, metrics_final.sHM_max, metrics_final.sHM_max/max(metrics_ref.sHM_max,eps), ...
        metrics_ref.u_max,   metrics_final.u_max,   metrics_final.u_max  /max(metrics_ref.u_max,  eps), ...
        volFrac);
end
