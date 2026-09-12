function saveStressSIMPHistoryCsv(resultRoot, history, configs)
    nRows = numel(history.iteration);
    if isfield(history, 'configNames') && size(history.stressAggregate, 2) == numel(history.configNames)
        constraintNames = cellstr(history.configNames(:));
    else
        constraintNames = cellfun(@(cfg) cfg.name, configs, 'UniformOutput', false);
    end

    rows = repmat(struct(), nRows, 1);
    for i = 1:nRows
        rows(i).iteration = history.iteration(i);
        rows(i).volumeFraction = history.volumeFraction(i);
        rows(i).change = history.change(i);
        rows(i).iterationTimeSec = history.iterationTimeSec(i);
        for k = 1:numel(constraintNames)
            name = matlab.lang.makeValidName(constraintNames{k});
            rows(i).(['S_' name]) = history.stressAggregate(i, k);
            rows(i).(['gStress_' name]) = history.constraint(i, k);
            rows(i).(['maxStress_' name]) = history.maxStress(i, k);
        end
    end
    writetable(struct2table(rows), fullfile(resultRoot, 'history.csv'));
end
