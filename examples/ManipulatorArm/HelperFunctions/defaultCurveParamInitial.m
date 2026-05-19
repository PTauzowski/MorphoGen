function params = defaultCurveParamInitial(arm)
% defaultCurveParamInitial  Torsion-inspired initial curve parameters.
%
%   anglePlusDeg / angleMinusDeg allow the two helix families to have
%   independent inclination angles.  Both are initialised to 45 degrees
%   (equal-angle symmetric design).  The legacy angleDeg field is kept so
%   that old saved param structs loaded from disk remain usable.

    params = struct();
    params.phaseFrac = 0.0;
    params.spacingFactor = 18.0;
    params.widthFactor = 0.9;
    params.helixPlusWeight = 0.8;
    params.anglePlusDeg = 45.0;
    params.helixMinusWeight = 0.8;
    params.angleMinusDeg = 45.0;
    params.axialWeight = 0.35;
    params.bendingWeight = 0.35;
    params.ringWeight = 0.2;
    params.ringSpacingFactor = 12.0;
    params.bevelRingWeight = 0.75;        % inclined joint face (stress concentration end)
    params.bevelRingWidthFactor = 3.0;   % wider than helix curves by default
    params.flatRingWeight = 0.5;          % perpendicular z=0 face
    params.flatRingWidthFactor = 3.0;
    params.middleRingWeight = 0.0;
    params.baseDensity = max(arm.mma.xminValue, 0.04);

    % Legacy fields — not in the optimisation vector but recognised as
    % fallbacks by buildCurveLinkedDensity when loading old result files.
    params.angleDeg = 45.0;
    params.jointRingWeight = 0.0;
    params.jointRingWidthFactor = params.widthFactor;
end
