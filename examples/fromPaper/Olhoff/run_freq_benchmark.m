% Benchmark runner for Olhoff & Du (2014) 2D beam cases.
% Runs CC, CS, SS boundary conditions and reports initial/final eigenfrequencies.
%
% Paper targets (Section 3.1):
%   CC: initial 146.1 -> optimal 456.4 (212% improvement)
%   CS: initial 104.1 -> optimal 288.7 (177% improvement)
%   SS: initial 68.7  -> optimal 174.7 (154% improvement)

clear; clc; close all;

L = 8; H = 1;
nelx = 240; nely = 30;
volFrac = 0.5;
penal = 3.0;
rmin  = 2*L/nelx;
maxiter = 300;    % enough for time-based beta continuation to reach beta=32
J = 3;

% Paper reference values
paper = struct();
paper.CC = struct('init', 146.1, 'opt', 456.4);
paper.CS = struct('init', 104.1, 'opt', 288.7);
paper.SS = struct('init', 68.7,  'opt', 174.7);

cases = { ...
    struct('code',"CC", 'label','Clamped–Clamped'); ...
    struct('code',"CS", 'label','Clamped–Simply'); ...
    struct('code',"SS", 'label','Simply–Simply'); ...
};

results = cell(numel(cases), 1);
for k = 1:numel(cases)
    c = cases{k};
    fprintf('\n================== %s ==================\n', c.label);
    fprintf('Paper: init=%.1f, opt=%.1f\n', paper.(c.code).init, paper.(c.code).opt);
    fprintf('=========================================\n');

    opts = struct('doDiagnostic',true,'diagnosticOnly',false,'diagModes',5);
    tic;
    [omega_best, xPhys_best, diag_out] = topFreqOptimization_MMA( ...
        L, H, nelx, nely, volFrac, penal, rmin, maxiter, c.code, J, opts);
    elapsed = toc;

    results{k} = struct('code', c.code, 'label', c.label, ...
        'diag', diag_out, 'omega', omega_best, 'xPhys', xPhys_best, 'time', elapsed);

    fprintf('\n--- %s Summary ---\n', c.label);
    print_block('Initial', diag_out.initial);
    print_block('Final',   diag_out.final);
    fprintf('Volume: %.4f (target 0.5)\n', mean(xPhys_best));
    fprintf('Grayness: %.4f\n', mean(4*xPhys_best.*(1-xPhys_best)));
    fprintf('Time: %.1f sec\n', elapsed);

    % Plot topology
    figure('Position', [100+300*(k-1), 100, 400, 100]);
    imagesc(1 - reshape(xPhys_best, nely, nelx));
    axis equal tight off; colormap(gray(256));
    title(sprintf('%s: ω₁=%.1f (paper: %.1f)', c.code, omega_best, paper.(c.code).opt));
end

% Final summary table
fprintf('\n\n========== FINAL BENCHMARK SUMMARY ==========\n');
fprintf('BC   | Init(code) | Init(paper) | Opt(code) | Opt(paper) | Improve\n');
fprintf('-----|------------|-------------|-----------|------------|--------\n');
for k = 1:numel(cases)
    r = results{k};
    p = paper.(r.code);
    improv = (r.diag.final.omega(1) / r.diag.initial.omega(1) - 1) * 100;
    fprintf('%s  | %10.1f | %11.1f | %9.1f | %10.1f | %5.0f%%\n', ...
        r.code, r.diag.initial.omega(1), p.init, r.diag.final.omega(1), p.opt, improv);
end
fprintf('==============================================\n');

function print_block(name, data)
    fprintf('%s eigenfreqs (rad/s): %8.2f %8.2f %8.2f\n', name, data.omega(1:min(3,end)));
end
