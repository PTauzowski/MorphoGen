function z_sharp = heavisideSharpenPost(z, Wfilter, beta, eta, fixedVars)
% HEAVISIDESHARPENPOST  Post-optimisation Heaviside sharpening of a density field.
%
%   z_sharp = heavisideSharpenPost(z, Wfilter, beta, eta, fixedVars)
%
%   Applies the density filter once, then projects through the smooth
%   Heaviside.  When beta is a vector, values are evaluated independently
%   and the result for the largest beta is returned.  The result is nearly
%   binary at large beta (>= 32) and is suitable for thresholding at 0.5.
%
%   This is a post-processing operation only — no re-optimisation is
%   performed.  For solveSIMPComplianceVolumeMMA the filter was applied
%   to the sensitivity during optimisation; applying it here as a density
%   filter is a valid approximation for boundary sharpening purposes.
%
%   Inputs:
%     z         [H x 1]  design variable (reference half-segment density)
%     Wfilter   [H x H]  sensitivity/density filter matrix
%     beta      scalar or [k x 1]  Heaviside sharpness parameter(s).
%               Values are evaluated independently; the largest beta is
%               returned.  Typical scalar calls: 8, 16, 32, 64.
%     eta       scalar  projection threshold in (0,1), default 0.5
%     fixedVars [m x 1]  indices forced to 1.0 after each projection step,
%               e.g. const ring elements.  Default: [] (none).
%
%   Output:
%     z_sharp   [H x 1]  sharpened density at the final (largest) beta

    if nargin < 4 || isempty(eta),       eta       = 0.5; end
    if nargin < 5 || isempty(fixedVars), fixedVars = [];  end

    beta = sort(beta(:), 'ascend');

    % Apply density filter once, then clamp to [0,1] (matches physicalDesign)
    z_filt = Wfilter * z(:);
    z_filt = min(max(z_filt, 0.0), 1.0);

    % Independent Heaviside projection; if beta is a vector, the sorted loop
    % leaves z_sharp at the largest beta.
    z_sharp = z_filt;
    for k = 1:numel(beta)
        z_sharp = applyHeaviside(z_filt, beta(k), eta);
        if ~isempty(fixedVars)
            z_sharp(fixedVars) = 1.0;
        end
    end
end

% -------------------------------------------------------------------------
function rho = applyHeaviside(rhoTilde, beta, eta)
    beta  = max(beta, eps);
    denom = tanh(beta * eta) + tanh(beta * (1.0 - eta));
    rho   = (tanh(beta * eta) + tanh(beta * (rhoTilde - eta))) / denom;
end
