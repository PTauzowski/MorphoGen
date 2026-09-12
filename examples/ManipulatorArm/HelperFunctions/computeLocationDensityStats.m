function stats = computeLocationDensityStats(model, x)
    nHalf = model.halfSegmentNelems;
    nLocations = floor(numel(x) / nHalf);
    assert(nLocations * nHalf == numel(x), ...
        'Full-arm element count %d is not divisible by half-segment count %d.', ...
        numel(x), nHalf);

    stats = repmat(struct('locationIndex', 0, 'elemStart', 0, 'elemEnd', 0, ...
        'meanDensity', 0, 'stdDensity', 0, 'rhoGt05', 0, 'rhoGt07', 0), ...
        nLocations, 1);
    for i = 1:nLocations
        ids = ((i - 1) * nHalf + 1):(i * nHalf);
        xi = x(ids);
        stats(i).locationIndex = i;
        stats(i).elemStart = ids(1);
        stats(i).elemEnd = ids(end);
        stats(i).meanDensity = mean(xi);
        stats(i).stdDensity = std(xi);
        stats(i).rhoGt05 = mean(xi > 0.5);
        stats(i).rhoGt07 = mean(xi > 0.7);
    end
end
