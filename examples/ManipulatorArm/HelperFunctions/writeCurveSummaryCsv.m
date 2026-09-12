function writeCurveSummaryCsv(resultRoot, params, bestMetrics, baselineMetrics, configs)
% writeCurveSummaryCsv  Save parameter and per-configuration metric summaries.

    names = curveParamNames();
    pRows = repmat(struct(), numel(names), 1);
    for i = 1:numel(names)
        pRows(i).parameter = string(names{i});
        pRows(i).value = params.(names{i});
    end
    writetable(struct2table(pRows), fullfile(resultRoot, 'best_parameters.csv'));

    rows = repmat(struct(), numel(configs), 1);
    for k = 1:numel(configs)
        rows(k).configName = string(configs{k}.name);
        rows(k).configLabel = string(configs{k}.label);
        rows(k).baselineMaxHM = baselineMetrics.maxHM(k);
        rows(k).bestMaxHM = bestMetrics.maxHM(k);
        rows(k).baselineStressAggregate = baselineMetrics.stressAggregateByConfig(k);
        rows(k).bestStressAggregate = bestMetrics.stressAggregateByConfig(k);
        rows(k).baselineTipUz = baselineMetrics.tipUz(k);
        rows(k).bestTipUz = bestMetrics.tipUz(k);
        rows(k).baselineMaxDisp = baselineMetrics.maxDisp(k);
        rows(k).bestMaxDisp = bestMetrics.maxDisp(k);
    end
    writetable(struct2table(rows), fullfile(resultRoot, 'metrics_by_config.csv'));
end
