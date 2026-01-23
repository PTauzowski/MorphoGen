%% ========================================================================
%  MMC Natural Frequency Maximization - Clamped/Clamped (Figure b)
%  ========================================================================
%
%  Design domain: 4 m x 1 m
%  Mesh         : 400 x 100
%  Supports     : Clamped on both ends
%  Mass         : At mid-span, mid-height
%
% =========================================================================

clear; clc; close all;

config.domainWidth  = 4;
config.domainHeight = 1;
config.numElemX     = 400;
config.numElemY     = 100;
config.boundaryType = 'clamped-both';

% Mass at mid-span
config.massLocation = [0, 0];

mmc_nfreq_2d_common(config);

