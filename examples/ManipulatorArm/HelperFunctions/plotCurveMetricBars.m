function fig = plotCurveMetricBars(metricsTable, metricName, titleText)
% plotCurveMetricBars  Compare variants by configuration for one metric.

    variants = unique(metricsTable.variant, 'stable');
    configs = unique(metricsTable.configLabel, 'stable');
    data = nan(numel(configs), numel(variants));
    for i = 1:numel(configs)
        for j = 1:numel(variants)
            row = metricsTable.configLabel == configs(i) & metricsTable.variant == variants(j);
            if any(row)
                data(i, j) = metricsTable.(metricName)(find(row, 1));
            end
        end
    end

    fig = figure('Color', 'white', 'Visible', 'on');
    b = bar(data);
    grid on;
    set(gca, 'XTick', 1:numel(configs), 'XTickLabel', cellstr(configs), ...
        'XTickLabelRotation', 20);
    legend(b(1:min(numel(b), numel(variants))), cellstr(variants(1:min(numel(b), numel(variants)))), ...
        'Location', 'best', 'Interpreter', 'none');
    ylabel(metricName, 'Interpreter', 'none');
    title(titleText, 'Interpreter', 'none');
end
