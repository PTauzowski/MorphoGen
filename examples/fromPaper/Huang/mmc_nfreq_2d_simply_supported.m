%% ========================================================================
%  MMC Natural Frequency Maximization - Simply Supported (Figure c)
%  ========================================================================
clear; clc; close all;

cfg = struct();

% --- Optimization ---
cfg.volfrac         = 0.30;
cfg.superEllipsoidPower = 6;
cfg.ksParameter     = 100;
cfg.maxIterations   = 500;
cfg.convergenceTol  = 1e-5;
cfg.moveLimit       = 0.10;
cfg.eps_initial     = 0.02;
cfg.eps_final       = 0.30;
cfg.eps_midIter     = 30;
cfg.eps_rate        = 0.50;

% --- Material ---
cfg.E0      = 2e11;
cfg.nu      = 0.3;
cfg.rho0    = 7800;
cfg.rho_min = 1e-6;
cfg.t       = 0.01;
cfg.m_lumped = 1e5;

% --- Geometry / mesh ---
cfg.L    = 4;
cfg.H    = 1;
cfg.nelx = 400;
cfg.nely = 100;

% --- Boundary condition & loads ---
cfg.boundaryType = 'simply-supported-both';
cfg.massLocation = [0, 0];

% --- Eigen solver ---
cfg.targetMode      = 1;
cfg.numModesCompute = 6;

% --- MMC initialization ---
cfg.componentSpacingX = 0.25;
cfg.componentSpacingY = 0.25;
cfg.initialHalfLength = 0.4;
cfg.initialHalfWidth1 = 0.04;
cfg.initialHalfWidth2 = 0.04;
cfg.initialAngle      = pi/4;

% --- Misc ---
cfg.seed = [];
cfg.verbosity = 1;

mmc_nfreq_2d_common(cfg);
