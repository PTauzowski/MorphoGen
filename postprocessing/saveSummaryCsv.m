function saveSummaryCsv(resultRoot, fd, history, finalX, finalC, VolFrac, ...
        objectiveDecreased, objectiveWindowConverged, mmaAcceptable, volumeActive, ...
        topologyNonuniform, topologyDiffersByArmLocation, ...
        locationDensityRange, locationDensityStd)
    row.status = "ok";
    row.allConfigurationsSolved = true;
    row.nElements = numel(finalX);
    row.nIterations = numel(history.iteration) - 1;
    row.fdMaxRelativeError = fd.maxRelativeError;
    row.fdPass = fd.maxRelativeError < 1.0e-2;
    row.initialJ = history.J(1);
    row.finalJ = history.J(end);
    row.relativeJChange = (history.J(end) - history.J(1)) / max(abs(history.J(1)), eps);
    row.targetVolumeFraction = VolFrac;
    row.finalVolumeFraction = mean(finalX);
    row.volumeActive = volumeActive;
    row.finalChange = history.change(end);
    row.objectiveDecreased = objectiveDecreased;
    row.objectiveWindowConverged = objectiveWindowConverged;
    row.mmaAcceptable = mmaAcceptable;
    row.topologyNonuniform = topologyNonuniform;
    row.finalDensityStd = std(finalX);
    row.finalDensityRange = max(finalX) - min(finalX);
    row.rhoGt05 = mean(finalX > 0.5);
    row.rhoGt07 = mean(finalX > 0.7);
    row.finalCMaxBending = finalC(1);
    row.finalCMaxTorsion = finalC(2);
    row.finalCMaxShear = finalC(3);
    row.topologyDiffersByArmLocation = topologyDiffersByArmLocation;
    row.locationMeanDensityRange = locationDensityRange;
    row.locationMeanDensityStd = locationDensityStd;
    writetable(struct2table(row), fullfile(resultRoot, 'summary.csv'));
end
