function plotLocationDensity(locationStats, resultRoot)
    loc = [locationStats.locationIndex]';
    meanDensity = [locationStats.meanDensity]';
    rhoGt05 = [locationStats.rhoGt05]';

    fig = figure('Name', 'Test A location density');
    tiledlayout(2, 1);

    nexttile;
    bar(loc, meanDensity);
    grid on; xlabel('Half-segment location'); ylabel('Mean density');

    nexttile;
    bar(loc, rhoGt05);
    grid on; xlabel('Half-segment location'); ylabel('Fraction \rho > 0.5');

    saveas(fig, fullfile(resultRoot, 'location_density.png'));
    savefig(fig, fullfile(resultRoot, 'location_density.fig'));
    close(fig);
end
