function exportInitialCurvesSTL(varargin)
% exportInitialCurvesSTL
% Write initial-curve topology STL for one module at volume fraction 0.5.
% No optimization is run. Output goes to the codex_smoke result folder.
%
% Usage:
%   exportInitialCurvesSTL()
%   exportInitialCurvesSTL('spacingFactor', 12)

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

spacingFactor = [];
for k = 1:2:numel(varargin)
    if strcmp(varargin{k}, 'spacingFactor')
        spacingFactor = varargin{k+1};
    end
end

resultRoot = fullfile(scriptDir, 'results', ...
    'testCurveParamRobustBetaOptimization_codex_smoke');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

%% Model
arm   = armModelDefaults("thin");
penal = 3.0;

nominalBeta = zeros(1, 7);
nominalBeta(2)=30;
fprintf('Building reference model...\n');
modelRef = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
    arm.res, arm.res_th, arm.alpha, nominalBeta, arm.ShapeFn, ...
    false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
fprintf('  Half-segment elements: %d\n', modelRef.halfSegmentNelems);

constRef  = armConstRingElementIds(modelRef, arm, "linked");
constFull = armConstRingElementIds(modelRef, arm, "full");

%% Initial density at Vf = 0.5
p0 = defaultCurveParamInitial(arm);
p0.axialWeight = 0.2;
p0.anglePlusDeg = 30.0;
p0.angleMinusDeg = 30.0;

% p0.helixPlusWeight = 0.8;
% p0.helixMinusWeight = 0.8;
% p0.axialWeight = 0.35;
% p0.bendingWeight = 0.35;
% p0.ringWeight = 0.5;
% p0.bevelRingWeight = 0.75;
% p0.flatRingWeight = 0.5;
% p0.middleRingWeight = 0.5;

p0.helixPlusWeight  = 0.5;
p0.helixMinusWeight = 0.5;
p0.axialWeight      = 0.5;
p0.bendingWeight    = 0.5;
p0.ringWeight       = 0.0;
p0.bevelRingWeight  = 0.5;
p0.flatRingWeight   = 0.5;
p0.middleRingWeight = 0.0;
   
p0.spacingFactor = helicesToSpacingFactor(4, arm);

if ~isempty(spacingFactor)
    p0.spacingFactor = spacingFactor;
end
fprintf('  spacingFactor = %.1f  (~%.0f helices around circumference)\n', ...
    p0.spacingFactor, 2*pi*arm.R / (p0.spacingFactor * arm.nominalElementSize));

opts = struct();
opts.penal          = penal;
opts.constRefElems  = constRef;
opts.constFullElems = constFull;
opts.rhoMin         = arm.mma.xminValue;
opts.elemSize       = arm.nominalElementSize;
opts.VolFrac        = 0.5;

fprintf('Building initial curve density at Vf=0.50...\n');
[~, rhoFull, info, fields] = buildCurveLinkedDensity(p0, modelRef, opts);
fprintf('  Achieved Vf = %.4f\n', info.fullVolumeFraction);

%% Export one-module STL
stlOut = fullfile(resultRoot, 'initial_curves_module.stl');
exportDensitySTL(modelRef, rhoFull, stlOut, 0.5, 1);
fprintf('Written: %s\n', stlOut);

[originRef, originNames] = classifyCurveFamilyOrigin(fields.reference, 0.85);
originFull = expandReferenceLabelsToArm(modelRef, originRef);

smoothStlOut = fullfile(resultRoot, 'initial_curves_module_smooth.stl');
smoothOpts = struct();
smoothOpts.nodalSmoothingIters = 1;
smoothOpts.surfaceSmoothingIters = 12;
smoothOpts.labels = originFull;
smoothOpts.labelNames = originNames;
smoothOpts.coloredObjPath = fullfile(resultRoot, 'initial_curves_module_colored.obj');
smoothInfo = exportDensitySmoothSTL(modelRef, rhoFull, smoothStlOut, 0.5, 1, smoothOpts);
fprintf('Written: %s  (%d vertices, %d triangles)\n', ...
    smoothStlOut, smoothInfo.nVertices, smoothInfo.nFaces);

objOut = smoothOpts.coloredObjPath;
fprintf('Written: %s  (smooth colored OBJ/MTL)\n', objOut);

voxelObjOut = fullfile(resultRoot, 'initial_curves_module_voxel_colored.obj');
objInfo = exportDensityColoredOBJ(modelRef, rhoFull, originFull, originNames, voxelObjOut, 0.5, 1);
fprintf('Written: %s  (%d vertices, %d triangles, %d material groups; exact voxel labels)\n', ...
    voxelObjOut, objInfo.nVertices, objInfo.nFaces, objInfo.nGroups);
end
