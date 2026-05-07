function bounds = defaultCurveParamBounds(arm)
% defaultCurveParamBounds  Bounds for curve-parametrized segment designs.
%
%   anglePlusDeg / angleMinusDeg are the per-family helix inclination bounds,
%   replacing the old shared angleDeg.  The legacy angleDeg entry is kept so
%   that old code calling curveParamsToUnconstrained on saved param structs
%   does not error if it tries to look up that field.

    elemSize = arm.nominalElementSize;
    bounds = struct();
    bounds.phaseFrac = [0.0, 1.0];               % circumferential phase as fraction of spacing
    bounds.spacingFactor = [8.0, 40.0];          % curve spacing / nominal element size
    bounds.widthFactor = [0.35, 2.2];            % curve Gaussian width / nominal element size
    bounds.helixPlusWeight = [0.0, 1.5];         % + helix family weight
    bounds.anglePlusDeg = [20.0, 70.0];          % + helix inclination from circumferential
    bounds.helixMinusWeight = [0.0, 1.5];        % - helix family weight
    bounds.angleMinusDeg = [20.0, 70.0];         % - helix inclination from circumferential
    bounds.axialWeight = [0.0, 0.8];             % straight axial family contribution
    bounds.bendingWeight = [0.0, 0.8];           % bending-plane side rib contribution
    bounds.ringWeight = [0.0, 0.4];              % periodic circumferential/ring contribution
    bounds.ringSpacingFactor = [12.0, 80.0];     % ring spacing / nominal element size
    bounds.jointRingWeight = [0.0, 1.5];         % extra reinforcement near segment ends
    bounds.jointRingWidthFactor = [0.5, 12.0];   % joint-ring Gaussian width / nominal element size
    bounds.baseDensity = [arm.mma.xminValue, min(0.35, 0.9 * arm.mma.xminValue + 0.30)];
    bounds.elemSize = [elemSize, elemSize];

    % Legacy — present for backward compatibility only, not in curveParamNames().
    bounds.angleDeg = [20.0, 70.0];
end
