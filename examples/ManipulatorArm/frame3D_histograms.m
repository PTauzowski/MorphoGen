%% coupledExamples.m
% Generate illustrative comparison figures for the coupled frame-solid
% solver and the corresponding full 3D solid analysis. The script exports
% paired deformed-mesh and Huber-Mises stress figures for representative
% joint-angle configurations and reports the frame cross-section constants
% used by the beam surrogate.

baseDir = fileparts(mfilename('fullpath'));
repoDir = fullfile(baseDir, '..', '..');

addpath(genpath(fullfile(repoDir, 'analysis')));
addpath(genpath(fullfile(repoDir, 'design')));
addpath(genpath(fullfile(repoDir, 'elements')));
addpath(genpath(fullfile(repoDir, 'examples', 'models')));
addpath(genpath(fullfile(repoDir, 'materials')));
addpath(genpath(fullfile(repoDir, 'math')));
addpath(genpath(fullfile(repoDir, 'mesh')));
addpath(genpath(fullfile(repoDir, 'examples', 'ManipulatorArm', 'HelperFunctions')));
addpath(genpath(fullfile(repoDir, 'postprocessing')));

docsDir = fullfile(repoDir, 'docs');
if ~exist(docsDir, 'dir')
    mkdir(docsDir);
end

close all;

arm = armModelDefaults("thin");
E = arm.E;
nu = arm.nu;
ShapeFn = arm.ShapeFn;
R = arm.R;
r = arm.r;
alpha = arm.alpha;
segmentLength = arm.segmentLength;
res = arm.res;
res_thickness = arm.res_th;
Pz = arm.Pz;

examples = struct( ...
    'name', { ...
        'sampleMinN_smooth', ...
        'sampleMaxMz_smooth', ...
        'sampleMaxTy_smooth', ...
        'sampleMaxMs_smooth' ...
    }, ...
    'sample', { ...
        [0 180 180 180 180 180 180], ...
        [0 0 0 180 180 180 180], ...
        [0 0 180 0 180 180 180], ...
        [0 45 45 45 270 180 180] ...
    }, ...
    'caption', { ...
        'Minimum axial-force configuration', ...
        'Maximum M_z configuration', ...
        'Maximum T_y configuration', ...
        'Maximum torsion M_s configuration' ...
    } ...
);

frameElems=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
nDiv  = 8;                           % 8 bins per angle
%vals  = linspace(0, 360 - 360/nDiv, nDiv);  % [0,45,90,...,315]
vals  = linspace(-180, 180 - 360/nDiv, nDiv);  % [0,45,90,...,315]

% We keep the first column fixed at 0, vary the remaining 6
nVar = nArms - 1;                    % 6
G = cell(1, nVar);
[G{:}] = ndgrid(vals);

% Assemble samples: rows = 8^6, cols = 7 (first col fixed to 0)
samples = zeros(nDiv^nVar, nArms);
%samples(:,1) = 0;
for k = 1:nVar
    samples(:, k+1) = G{k}(:);
end
nSamples = size(samples,1);

% same sample count, continuous uniform draw
samples(:, 2:end) = -180 + 360 * lhsdesign(nSamples, nVar);


%[~, ~, frameNodes] = computeArmSamples(E, nu, segmentLength, alpha, samples );
%[vN, vTy, vTz, vMs, vMy, vMz, all_forces] = computeAllInternalForces( frameNodes, frameElems, E, nu, R, r, samples, Pz);

%save("ManipulatorFrame3Ddata200K.mat");
load("ManipulatorFrame3Ddata200K.mat");

nbin = 2000;

figsDir = fullfile(baseDir, 'figs');
plotInternalForceHistograms(vN, vTy, vTz, vMs, vMy, vMz, figsDir, nbin, false);
plotSnakeView(vN, vTy, vTz, vMs, vMy, vMz, figsDir, false);
plotSnakeViewNodes(vN, vTy, vTz, vMs, vMy, vMz, figsDir, false);
