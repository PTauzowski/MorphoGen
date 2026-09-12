function sf = helicesToSpacingFactor(nHelices, arm)
% helicesToSpacingFactor  Convert number of helices around circumference to spacingFactor.
%
%   sf = helicesToSpacingFactor(nHelices, arm)
%
%   Inverse of the display formula used throughout the arm pipeline:
%     nHelices = 2*pi*arm.R / (spacingFactor * arm.nominalElementSize)
%
%   nHelices - desired number of helix curves around the circumference
%   arm      - arm model struct from armModelDefaults (needs R, nominalElementSize)

sf = 2*pi * arm.R / (nHelices * arm.nominalElementSize);
end
