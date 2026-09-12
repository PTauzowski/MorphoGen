function plotInternalForceHistograms(vN, vTy, vTz, vMs, vMy, vMz, figsDir, nbin, skipNoise)
% plotInternalForceHistograms  Save histograms (up to 12 forces × 6 segments).
%
%   plotInternalForceHistograms(vN,vTy,vTz,vMs,vMy,vMz,figsDir,nbin,skipNoise)
%
%   Each force array is (nSamples x 2 x 6):
%     dimension 2 : 1 = start node, 2 = end node
%     dimension 3 : segment index 1-6
%
%   skipNoise (default true):
%     true  — channels where max(|data|) < NOISE_THRESH are skipped
%     false — all 72 histograms are saved regardless
%
%   PNG files are written to figsDir with names like:
%     hist_seg1_node1_N.png
%     hist_seg3_node2_Mz.png  ...

% Anything below this is indistinguishable from FEM floating-point noise.
% Physical forces here are O(100 N) / O(100 N·m); noise is O(1e-10).
NOISE_THRESH = 1e-6;

if nargin < 8 || isempty(nbin)
    nbin = 500;
end
if nargin < 9 || isempty(skipNoise)
    skipNoise = true;
end

if ~exist(figsDir, 'dir')
    mkdir(figsDir);
end

forces = { vN,  vTy,  vTz,  vMs,  vMy,  vMz  };
fkeys  = { 'N', 'Ty', 'Tz', 'Ms', 'My', 'Mz' };
flabels = { ...
    'Axial force N (N)', ...
    'Shear T_y (N)', ...
    'Shear T_z (N)', ...
    'Torsional moment M_s (N·m)', ...
    'Bending moment M_y (N·m)', ...
    'Bending moment M_z (N·m)'  };

nodeLabels = {'start node (1)', 'end node (2)'};

nSegs   = 6;
nNodes  = 2;
nSkip   = 0;
nSaved  = 0;

fig = figure('Visible', 'off');

for seg = 1:nSegs
    for nd = 1:nNodes
        for fi = 1:numel(forces)
            data = forces{fi}(:, nd, seg);

            if skipNoise && max(abs(data)) < NOISE_THRESH
                fprintf('  Skipped (numerical noise): seg%d node%d %s  [max=%.2e]\n', ...
                    seg, nd, fkeys{fi}, max(abs(data)));
                nSkip = nSkip + 1;
                continue;
            end

            clf(fig);
            histogram(data, nbin);

            titleStr = sprintf('Segment %d — %s — %s', ...
                seg, flabels{fi}, nodeLabels{nd});
            title(titleStr, 'Interpreter', 'tex');
            xlabel(flabels{fi}, 'Interpreter', 'tex');
            ylabel('Count');
            grid on;

            fname = sprintf('hist_seg%d_node%d_%s', seg, nd, fkeys{fi});
            exportgraphics(fig, fullfile(figsDir, [fname '.png']), 'Resolution', 150);
            set(fig, 'Visible', 'on');
            savefig(fig, fullfile(figsDir, [fname '.fig']));
            set(fig, 'Visible', 'off');
            nSaved = nSaved + 1;
        end
    end
end

close(fig);
fprintf('Histograms saved: %d   Skipped (noise): %d\n', nSaved, nSkip);
end
