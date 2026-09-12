function saveStressSIMPSummaryCsv(resultRoot, history, finalX, finalStressAggregate, ...
        finalMaxStress, finalConstraint, stressCoeff, stressPNorm, stressRelaxationQ, ...
        constraintsSatisfied, topologyNonuniform, topologyDiffersByArmLocation, ...
        locationDensityRange, locationDensityStd)
    row.status = "ok";
    row.nElements = numel(finalX);
    row.nIterations = numel(history.iteration) - 1;
    row.stressCoeff = stressCoeff;
    row.stressPNorm = stressPNorm;
    row.stressRelaxationQ = stressRelaxationQ;
    row.finalVolumeFraction = mean(finalX);
    row.finalChange = history.change(end);
    row.maxStressConstraint = max(finalConstraint);
    row.constraintsSatisfied = constraintsSatisfied;
    row.topologyNonuniform = topologyNonuniform;
    row.finalDensityStd = std(finalX);
    row.finalDensityRange = max(finalX) - min(finalX);
    row.rhoGt05 = mean(finalX > 0.5);
    row.rhoGt07 = mean(finalX > 0.7);
    row.finalStressAggregateMax = max(finalStressAggregate);
    row.finalMaxStress = max(finalMaxStress);
    row.topologyDiffersByArmLocation = topologyDiffersByArmLocation;
    row.locationMeanDensityRange = locationDensityRange;
    row.locationMeanDensityStd = locationDensityStd;
    writetable(struct2table(row), fullfile(resultRoot, 'summary.csv'));
end
