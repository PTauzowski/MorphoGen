function arm = armModelDefaults(preset)
% armModelDefaults  Shared material, geometry, and optimizer defaults.
%
%   arm = armModelDefaults("thick") returns the legacy thick-wall model.
%
%   arm = armModelDefaults("thin") returns the thin-wall publication model.
    if nargin < 1
        preset = "thin";
    end
    preset = string(validatestring(char(preset), {'thick', 'thin'}));

    arm = struct();
    arm.preset = preset;

    arm.E = 2.0e9;
    arm.nu = 0.35;
    arm.R = 0.14;
    arm.alpha = 22.5;
    arm.alpha_deg = arm.alpha;
    arm.res = 15;
    arm.Pz = 100;
    arm.ShapeFn = ShapeFunctionH8();

    switch preset
        case "thick"
            arm.r = 0.08;
            arm.h_seg = 0.25;
            arm.res_th = 4;
        case "thin"
            arm.r = 0.13;
            arm.h_seg = 0.25;
            arm.res_th = 1;
    end

    arm.h = arm.h_seg;
    arm.segmentLength = arm.h_seg;
    arm.wallThickness = arm.R - arm.r;
    arm.nominalElementSize = arm.wallThickness / arm.res_th;
    arm.Rfilter = 3 * arm.nominalElementSize;
    arm.Rmin = arm.Rfilter;
    arm.filletR = 0;   % joint fillet radius [m]; 0 = sharp (no fillet)

    % Constant-density ring controls. Each enabled ring is one FE element
    % wide along the local half-segment axis.
    arm.constEndRing = true;
    arm.constMiddleRing = false;

    % Circumferential divisor for rotation-aware linking (Option 1).
    % resCirc is snapped to the nearest multiple of nCircDiv so that every
    % joint angle used in armLoadConfigs("six") is an exact integer multiple
    % of the mesh angular step 2*pi/resCirc.
    % nCircDiv = 8 covers all multiples of 45° (GCD of {45,90,180,270}).
    arm.nCircDiv = 8;

    arm.mma = struct();
    arm.mma.xminValue = 0.01;
    arm.mma.changeTol = 1.0e-3;
    arm.mma.moveLimit = 0.05;
    arm.mma.minMoveLimit = 0.003;
    arm.mma.moveDecay = 0.95;
    arm.mma.mmaDamping = 0.25;
end
