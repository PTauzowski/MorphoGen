function plotSnakeViewNodes(vN, vTy, vTz, vMs, vMy, vMz, figsDir, skipNoise)
% plotSnakeViewNodes  Node-centric snake view showing discontinuity at joints.
%
%   plotSnakeViewNodes(vN,vTy,vTz,vMs,vMy,vMz,figsDir,skipNoise)
%
%   Each force array is (nSamples x 2 x nElems):
%     dimension 2 : 1 = start node, 2 = end node
%     dimension 3 : element index
%
%   X-axis runs over nodes (base=1 to tip=nElems+1).
%   At each internal node two box-whisker distributions are drawn:
%     orange = start node of outgoing element (local frame of elem i+1)
%     blue   = end node of incoming element   (local frame of elem i)
%   A dotted connector between them shows the magnitude of the discontinuity.
%   Terminal nodes (base and tip) carry only one distribution.
%
%   skipNoise (default true): skip channels where max(|data|) < 1e-6.
%
%   Produces one PNG per force type (6 total):
%     nodeSnake_N.png, nodeSnake_My.png, ...

NOISE_THRESH = 1e-6;
DELTA        = 0.18;   % x-offset for paired boxes at internal nodes
HW           = 0.12;   % half-width of each box

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

nElems = size(vN, 3);
nNodes = nElems + 1;

pcts = [5 25 50 75 95];

colOut = [0.84 0.33 0.10];   % outgoing: start node of element (orange)
colIn  = [0.18 0.44 0.79];   % incoming: end node of element   (blue)
colGap = [0.50 0.50 0.50];   % dotted discontinuity connector  (grey)

nSaved = 0;
nSkip  = 0;

fig = figure('Visible', 'off', 'Position', [100 100 1000 440]);

for fi = 1:numel(forces)

    % --- noise check ---
    allNoise = true;
    for seg = 1:nElems
        for nd = 1:2
            if max(abs(forces{fi}(:, nd, seg))) >= NOISE_THRESH
                allNoise = false; break;
            end
        end
        if ~allNoise, break; end
    end
    if skipNoise && allNoise
        fprintf('  Skipped (noise): nodeSnake %s\n', fkeys{fi});
        nSkip = nSkip + 1;
        continue;
    end

    % --- percentiles: Pout(p, elem) = start nodes, Pin(p, elem) = end nodes ---
    Pout = zeros(numel(pcts), nElems);
    Pin  = zeros(numel(pcts), nElems);
    for seg = 1:nElems
        Pout(:, seg) = prctile(forces{fi}(:, 1, seg), pcts);
        Pin(:, seg)  = prctile(forces{fi}(:, 2, seg), pcts);
    end

    % --- x positions ---
    % Outgoing (start nodes): elem 1 at node 1 (no offset), elems 2..nElems offset right
    xOut = [1,  (2:nElems) + DELTA];
    % Incoming (end nodes): elems 1..nElems-1 offset left, last elem at node nNodes (no offset)
    xIn  = [(2:nElems) - DELTA,  nNodes];

    clf(fig);
    hold on;

    % --- whisker boxes ---
    for seg = 1:nElems
        whiskerBox(xOut(seg), Pout(:,seg), colOut, HW);
        whiskerBox(xIn(seg),  Pin(:,seg),  colIn,  HW);
    end

    % --- median connection lines ---
    plot(xOut, Pout(3,:), '-o', 'Color', colOut, 'LineWidth', 1.8, ...
        'MarkerFaceColor', colOut, 'MarkerSize', 5, ...
        'DisplayName', 'Elem start node (outgoing)');
    plot(xIn,  Pin(3,:),  '-s', 'Color', colIn,  'LineWidth', 1.8, ...
        'MarkerFaceColor', colIn,  'MarkerSize', 5, ...
        'DisplayName', 'Elem end node (incoming)');

    % --- dotted discontinuity connectors at internal nodes ---
    for seg = 1:nElems-1
        xL = xIn(seg);          % end of elem seg   (left side of joint)
        xR = xOut(seg+1);       % start of elem seg+1 (right side of joint)
        mL = Pin(3, seg);
        mR = Pout(3, seg+1);
        plot([xL xR], [mL mR], ':', 'Color', colGap, ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
    end

    yline(0, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    hold off;

    xticks(1:nNodes);
    xticklabels(arrayfun(@(n) sprintf('N%d', n), 1:nNodes, 'UniformOutput', false));
    xlim([0.5,  nNodes + 0.5]);
    xlabel('Node  (base \rightarrow tip)');
    ylabel(flabels{fi}, 'Interpreter', 'tex');
    title(sprintf('Node snake view — %s  (local element coordinates)', flabels{fi}), ...
        'Interpreter', 'tex');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    box on;

    fname = sprintf('nodeSnake_%s', fkeys{fi});
    exportgraphics(fig, fullfile(figsDir, [fname '.png']), 'Resolution', 150);
    set(fig, 'Visible', 'on');
    savefig(fig, fullfile(figsDir, [fname '.fig']));
    set(fig, 'Visible', 'off');
    nSaved = nSaved + 1;
end

close(fig);
fprintf('Node snake plots saved: %d   Skipped (noise): %d\n', nSaved, nSkip);
end

% -------------------------------------------------------------------------
function whiskerBox(x, P, col, hw)
% Draw a single box-whisker: P = [p5 p25 median p75 p95], hw = half-width.
p5  = P(1);  p25 = P(2);  med = P(3);  p75 = P(4);  p95 = P(5);

% 5–95 whisker stem
plot([x x], [p5 p95], '-', 'Color', col, 'LineWidth', 0.8, 'HandleVisibility', 'off');
% Whisker caps
capW = hw * 0.45;
plot([x-capW x+capW], [p5  p5 ], '-', 'Color', col, 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot([x-capW x+capW], [p95 p95], '-', 'Color', col, 'LineWidth', 0.8, 'HandleVisibility', 'off');
% 25–75 filled box
patch([x-hw x+hw x+hw x-hw], [p25 p25 p75 p75], col, ...
    'FaceAlpha', 0.30, 'EdgeColor', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
% Median bar
plot([x-hw x+hw], [med med], '-', 'Color', col*0.55, ...
    'LineWidth', 2.2, 'HandleVisibility', 'off');
end
