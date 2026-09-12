function plotSnakeView(vN, vTy, vTz, vMs, vMy, vMz, figsDir, skipNoise)
% plotSnakeView  12 snake-view plots of internal forces along the arm.
%
%   plotSnakeView(vN,vTy,vTz,vMs,vMy,vMz,figsDir,skipNoise)
%
%   Each force array is (nSamples x 2 x 6):
%     dimension 2 : 1 = start node, 2 = end node
%     dimension 3 : segment index 1-6
%
%   skipNoise (default true):
%     true  — channels where all 6 segments are numerical noise are skipped
%     false — all 12 snake plots are saved regardless
%
%   One plot is produced per force type (N, Ty, Tz, Ms, My, Mz) per node
%   (start / end) = 12 plots total.  Each plot has 6 x-positions (one per
%   segment) and shows median + 25–75 pct band + 5–95 pct band.
%
%   PNG files written to figsDir:
%     snake_N_node1.png, snake_Mz_node2.png, ...

NOISE_THRESH = 1e-6;

if nargin < 8 || isempty(skipNoise)
    skipNoise = true;
end

if ~exist(figsDir, 'dir')
    mkdir(figsDir);
end

forces  = { vN,  vTy,  vTz,  vMs,  vMy,  vMz  };
fkeys   = { 'N', 'Ty', 'Tz', 'Ms', 'My', 'Mz' };
flabels = { ...
    'Axial force N (N)', ...
    'Shear T_y (N)', ...
    'Shear T_z (N)', ...
    'Torsional moment M_s (N{\cdot}m)', ...
    'Bending moment M_y (N{\cdot}m)', ...
    'Bending moment M_z (N{\cdot}m)' };

nodeLabels = {'start node (1)', 'end node (2)'};
nodeKeys   = {'node1', 'node2'};

nSegs = size(vN, 3);
xpos  = 1:nSegs;
xLabels = arrayfun(@(s) sprintf('Seg %d', s), xpos, 'UniformOutput', false);

pcts = [5 10 25 40 50 60 75 90 95];

colBand4 = [0.88 0.93 1.00];   % 5–95 pct fill
colBand3 = [0.67 0.80 0.96];   % 10–90 pct fill
colBand2 = [0.38 0.62 0.90];   % 25–75 pct fill
colBand1 = [0.15 0.38 0.75];   % 40–60 pct fill
colMed   = [0.05 0.18 0.52];   % median line

nSaved = 0;
nSkip  = 0;

fig = figure('Visible', 'off', 'Position', [100 100 860 420]);

for nd = 1:2
    for fi = 1:numel(forces)

        % Compute percentiles at each segment
        P = zeros(numel(pcts), nSegs);
        allNoise = true;
        for seg = 1:nSegs
            data = forces{fi}(:, nd, seg);
            if max(abs(data)) >= NOISE_THRESH
                allNoise = false;
            end
            P(:, seg) = prctile(data, pcts);
        end

        if skipNoise && allNoise
            fprintf('  Skipped (noise): snake %s %s\n', fkeys{fi}, nodeKeys{nd});
            nSkip = nSkip + 1;
            continue;
        end

        p5  = P(1,:);  p10 = P(2,:);  p25 = P(3,:);  p40 = P(4,:);
        med = P(5,:);
        p60 = P(6,:);  p75 = P(7,:);  p90 = P(8,:);  p95 = P(9,:);

        clf(fig);
        hold on;

        % 5–95 percentile band (outermost)
        fill([xpos, fliplr(xpos)], [p5, fliplr(p95)], colBand4, ...
            'EdgeColor', 'none', 'DisplayName', '5–95 pct');

        % 10–90 percentile band
        fill([xpos, fliplr(xpos)], [p10, fliplr(p90)], colBand3, ...
            'EdgeColor', 'none', 'DisplayName', '10–90 pct');

        % 25–75 percentile band
        fill([xpos, fliplr(xpos)], [p25, fliplr(p75)], colBand2, ...
            'EdgeColor', 'none', 'DisplayName', '25–75 pct');

        % 40–60 percentile band (innermost)
        fill([xpos, fliplr(xpos)], [p40, fliplr(p60)], colBand1, ...
            'EdgeColor', 'none', 'DisplayName', '40–60 pct');

        % Median
        plot(xpos, med, '-o', 'Color', colMed, 'LineWidth', 2.0, ...
            'MarkerFaceColor', colMed, 'MarkerSize', 6, 'DisplayName', 'Median');

        % Zero reference
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');

        hold off;

        xlim([0.6  nSegs + 0.4]);
        xticks(xpos);
        xticklabels(xLabels);
        xlabel('Segment (base \rightarrow tip)');
        ylabel(flabels{fi}, 'Interpreter', 'tex');
        title(sprintf('Snake view — %s — %s', flabels{fi}, nodeLabels{nd}), ...
            'Interpreter', 'tex');
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        box on;

        fname = sprintf('snake_%s_%s', fkeys{fi}, nodeKeys{nd});
        exportgraphics(fig, fullfile(figsDir, [fname '.png']), 'Resolution', 150);
        set(fig, 'Visible', 'on');
        savefig(fig, fullfile(figsDir, [fname '.fig']));
        set(fig, 'Visible', 'off');
        nSaved = nSaved + 1;
    end
end

close(fig);
fprintf('Snake plots saved: %d   Skipped (noise): %d\n', nSaved, nSkip);
end
