function saveHistoryCsv(resultRoot, history, configs)
    nRows = numel(history.iteration);
    rows = repmat(struct(), nRows, 1);
    for i = 1:nRows
        rows(i).iteration = history.iteration(i);
        rows(i).J = history.J(i);
        rows(i).volumeFraction = history.volumeFraction(i);
        rows(i).change = history.change(i);
        if isfield(history, 'iterationTimeSec')
            rows(i).iterationTimeSec = history.iterationTimeSec(i);
        end
        for k = 1:numel(configs)
            fieldName = ['C_' configs{k}.name];
            rows(i).(fieldName) = history.C(i, k);
        end
    end
    writetable(struct2table(rows), fullfile(resultRoot, 'history.csv'));
end
