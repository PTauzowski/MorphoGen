% Generate some results
x = linspace(0, 2*pi, 100);
y = sin(x);
z = cos(x);

% Create a figure for the plot
fig = figure;
plot(x, y, 'r', 'LineWidth', 1.5); hold on;
plot(x, z, 'b--', 'LineWidth', 1.5);
title('Sine and Cosine Functions');
xlabel('x');
ylabel('y');
legend('sin(x)', 'cos(x)');
grid on;

% Save the plot as a PDF
exportgraphics(fig, 'results_plot.pdf', 'ContentType', 'vector'); % Save plot only
close(fig);

% Add Text Results to PDF (Using Figures with Text)
result_fig = figure('Visible', 'off');
annotation('textbox', [0.1, 0.5, 0.8, 0.4], 'String', ...
    sprintf(['Computation Results:\n\n', ...
    'Max of sin(x): %.2f\n', ...
    'Min of sin(x): %.2f\n', ...
    'Max of cos(x): %.2f\n', ...
    'Min of cos(x): %.2f'], ...
    max(y), min(y), max(z), min(z)), ...
    'FontSize', 12, 'EdgeColor', 'none');
exportgraphics(result_fig, 'results_plot.pdf', 'ContentType', 'vector','Append', true);
close(result_fig);

% Combine PDFs using a command-line tool (if needed)
disp('PDFs created: results_plot.pdf and results_text.pdf');



% Define the PDF file name
pdfFileName = 'combined_plots.pdf';

% First Plot
fig1 = figure;
plot(1:10, rand(1, 10), '-o', 'LineWidth', 1.5);
title('Random Data 1');
xlabel('X');
ylabel('Y');
exportgraphics(fig1, pdfFileName, 'ContentType', 'vector', 'Append', false); % Create the PDF
close(fig1);

% Second Plot
fig2 = figure;
plot(1:10, rand(1, 10), '-x', 'LineWidth', 1.5);
title('Random Data 2');
xlabel('X');
ylabel('Y');
exportgraphics(fig2, pdfFileName, 'ContentType', 'vector', 'Append', true); % Append to the PDF
close(fig2);

% Third Plot
fig3 = figure;
plot(1:10, rand(1, 10), '-s', 'LineWidth', 1.5);
title('Random Data 3');
xlabel('X');
ylabel('Y');
exportgraphics(fig3, pdfFileName, 'ContentType', 'vector', 'Append', true); % Append to the PDF
close(fig3);

disp(['Combined PDF created: ', pdfFileName]);
