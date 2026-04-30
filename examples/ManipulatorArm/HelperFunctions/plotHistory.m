function plotHistory(history, configs, resultRoot)
    fig = figure('Name', 'Test A convergence');
    tiledlayout(3, 1);

    nexttile;
    plot(history.iteration, history.J, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('J');

    nexttile;
    plot(history.iteration, history.C, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('Compliance C_k');
    legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), ...
        'Location', 'best');

    nexttile;
    plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('Volume fraction');

    saveas(fig, fullfile(resultRoot, 'convergence_history.png'));
    close(fig);
end
