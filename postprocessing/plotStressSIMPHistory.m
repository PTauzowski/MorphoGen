function plotStressSIMPHistory(history, configs, resultRoot, testName)
    labels = cellfun(@(s) s.label, configs, 'UniformOutput', false);

    fig = figure('Name', testName + " volume history");
    plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.2);
    grid on;
    xlabel('Iteration');
    ylabel('Volume fraction');
    title(testName + ": volume minimization");
    saveas(fig, fullfile(resultRoot, 'volume_history.png'));
    savefig(fig, fullfile(resultRoot, 'volume_history.fig'));
    close(fig);

    fig = figure('Name', testName + " stress constraint history");
    plot(history.iteration, history.constraint, '-o', 'LineWidth', 1.0);
    yline(0, '--r', 'constraint limit');
    grid on;
    xlabel('Iteration');
    ylabel('Stress aggregate constraint g');
    legend(labels, 'Location', 'best');
    title(testName + ": stress constraints");
    saveas(fig, fullfile(resultRoot, 'stress_constraint_history.png'));
    savefig(fig, fullfile(resultRoot, 'stress_constraint_history.fig'));
    close(fig);

    fig = figure('Name', testName + " max stress history");
    plot(history.iteration, history.maxStress, '-o', 'LineWidth', 1.0);
    grid on;
    xlabel('Iteration');
    ylabel('Max HM stress');
    legend(labels, 'Location', 'best');
    title(testName + ": max stress");
    saveas(fig, fullfile(resultRoot, 'max_stress_history.png'));
    savefig(fig, fullfile(resultRoot, 'max_stress_history.fig'));
    close(fig);
end
