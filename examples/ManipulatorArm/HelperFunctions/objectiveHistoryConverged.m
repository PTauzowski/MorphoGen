function tf = objectiveHistoryConverged(JHistory, tol)
    objective = JHistory(~isnan(JHistory));
    nWindow = min(10, numel(objective));
    if nWindow < 3
        tf = false;
        return;
    end
    j = objective(end-nWindow+1:end);
    tf = (max(j) - min(j)) / max(abs(j(end)), eps) < tol;
end
