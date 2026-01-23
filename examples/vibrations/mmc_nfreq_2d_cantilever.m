%% ========================================================================
%  MMC Natural Frequency Maximization - Cantilever (Figure a)
%  ========================================================================
%
%  Design domain: 2 m x 1 m
%  Mesh         : 200 x 100
%  Supports     : Clamped on the left edge
%  Mass         : At the right tip, mid-height
%
% =========================================================================

clear; clc; close all;

config.domainWidth  = 2;
config.domainHeight = 1;
config.numElemX     = 200;
config.numElemY     = 100;
config.boundaryType = 'clamped-left';

% Mass at the tip (right edge, mid-height)
config.massLocation = [config.domainWidth/2, 0];

mmc_nfreq_2d_common(config);

