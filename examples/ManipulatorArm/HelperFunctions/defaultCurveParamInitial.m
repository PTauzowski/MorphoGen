function params = defaultCurveParamInitial(arm)
% defaultCurveParamInitial  Torsion-inspired initial curve parameters.

    params = struct();
    params.angleDeg = 45.0;
    params.phaseFrac = 0.0;
    params.spacingFactor = 18.0;
    params.widthFactor = 0.9;
    params.helixPlusWeight = 0.8;
    params.helixMinusWeight = 0.8;
    params.axialWeight = 0.35;
    params.bendingWeight = 0.35;
    params.ringWeight = 0.0;
    params.ringSpacingFactor = 32.0;
    params.jointRingWeight = 0.75;
    params.baseDensity = max(arm.mma.xminValue, 0.04);
end
