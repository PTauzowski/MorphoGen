% Benchmark runner for Olhoff & Du (2014) 2D beam cases.
% Runs CC, CS, SS boundary conditions and reports initial/final eigenfrequencies.

clear; clc; close all;

% Paper reference values
paper = struct();
paper.CC = struct('init', 146.1, 'opt', 456.4);
paper.CS = struct('init', 104.1, 'opt', 288.7);
paper.SS = struct('init', 68.7,  'opt', 174.7);


L = 8; H = 1;
nelx = 240; nely = 30;
volFrac = 0.5;
penal = 3.0;
rmin  = 2*L/nelx;
maxiter = 300;    % enough for time-based beta continuation to reach beta=32
J = 3;
 opts = struct('doDiagnostic',true,'diagnosticOnly',false,'diagModes',5);

% Simple supported beam
[omega_best, xPhys_best, diag_out] = topFreqOptimization_MMA( ...
        L, H, nelx, nely, volFrac, penal, rmin, maxiter, "SS", J, opts);
 title(sprintf('%s: ω₁=%.1f (paper: %.1f)', "SS", omega_best, paper.("SS").opt));

% Left side clamped  beam
figure, hold on;
[omega_best, xPhys_best, diag_out] = topFreqOptimization_MMA( ...
        L, H, nelx, nely, volFrac, penal, rmin, maxiter, "CS", J, opts);
 title(sprintf('%s: ω₁=%.1f (paper: %.1f)', "CS", omega_best, paper.("CS").opt));

% Clamped-clamped beam.
figure, hold on;
[omega_best, xPhys_best, diag_out] = topFreqOptimization_MMA( ...
        L, H, nelx, nely, volFrac, penal, rmin, maxiter, "CC", J, opts);
 title(sprintf('%s: ω₁=%.1f (paper: %.1f)', "CC", omega_best, paper.("CC").opt));