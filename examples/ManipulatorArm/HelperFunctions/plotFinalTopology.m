function plotFinalTopology(model, x, resultRoot)
    fig = figure('Name', 'Test A final density');
    plotElementDensityField(model, x);
    title(sprintf('Test A final density, vf=%.3f', mean(x)));
    saveas(fig, fullfile(resultRoot, 'final_density.png'));
    close(fig);

    thresholds = [0.5, 0.7];
    for i = 1:numel(thresholds)
        t = thresholds(i);
        stem = sprintf('final_threshold_rho_gt_%02d', round(10 * t));
        fig = figure('Name', sprintf('Test A rho > %.1f', t));
        hold on; axis on; daspect([1 1 1]); view(45, 35);
        selected = x > t;
        model.fe.plotSolidSelected(model.mesh.nodes, selected, [0.15 0.15 0.15]);
        xlabel('x'); ylabel('y'); zlabel('z');
        title(sprintf('Test A final topology, \\rho > %.1f', t));
        saveas(fig, fullfile(resultRoot, [stem '.png']));
        close(fig);
    end
end
