%% ========================================================================
%  MMC Natural Frequency Maximization - Cantilever (Figure a)
%  ========================================================================
clear; clc; close all;

cfg = struct();

% --- Optimization ---
cfg.volfrac         = 0.40;
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
cfg.L    = 2;
cfg.H    = 1;
cfg.nelx = 200;
cfg.nely = 100;

% --- Boundary condition & loads ---
cfg.boundaryType = 'clamped-left';
cfg.massLocation = [cfg.L/2, 0];   % tip mass at mid-height

% --- Eigen solver ---
cfg.targetMode      = 1;
cfg.numModesCompute = 6;

% --- MMC initialization (default grid kept explicit for reproducibility) ---
cfg.componentSpacingX = 0.25;
cfg.componentSpacingY = 0.25;
cfg.initialHalfLength = 0.4;
cfg.initialHalfWidth1 = 0.04;
cfg.initialHalfWidth2 = 0.04;
cfg.initialAngle      = pi/4;

% --- Misc ---
cfg.seed = [];      % set numeric value for deterministic runs
cfg.verbosity = 1;

mmc_nfreq_2d_common(cfg);
