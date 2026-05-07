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
    params.ringWeight = 0.0;
    params.ringSpacingFactor = 32.0;
    params.jointRingWeight = 0.75;
    params.jointRingWidthFactor = 3.0;   % wider than helix curves by default
    params.baseDensity = max(arm.mma.xminValue, 0.04);

    % Legacy field — not in the optimisation vector but recognised as
    % fallback by buildCurveLinkedDensity when loading old result files.
    params.angleDeg = 45.0;
end
