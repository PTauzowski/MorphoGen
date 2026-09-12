function fd = finiteDifferenceGradientTestLinked(model, analyses, rho, penal, pAgg, weights, C0, dJdrho, nTest, h, xmin, xmax, useParallel)
% FINITEDIFFERENCEGRADIENTTESTLINKED  Finite-difference gradient check in rho-space.
%
%   Validates the pullback gradient dJdrho by perturbing individual rho(e)
%   entries, expanding the full-arm density via model.segmentToArm, and
%   comparing the central-difference estimate of dJ/drho(e) against the
%   analytic value.
%
%   fd = finiteDifferenceGradientTestLinked(model, analyses, rho, penal,
%           pAgg, weights, C0, dJdrho, nTest, h, xmin, xmax)
%
%   Inputs:
%     model    ManipulatorModel3D  used for segmentToArm expansion
%     analyses {nConfigs x 1}     FE analysis objects for all configurations
%     rho      [H x 1]            current reference half-segment densities
%     penal    scalar             SIMP penalization exponent
%     pAgg     scalar             p-norm aggregation exponent
%     weights  [nConfigs x 1]    per-configuration weights
%     C0       [nConfigs x 1]    compliance normalization values
%     dJdrho   [H x 1]           analytic gradient in rho-space (from pullback)
%     nTest    scalar             number of rho indices to test
%     h        scalar             nominal finite-difference step size
%     xmin     [H x 1]           lower bounds on rho
%     xmax     [H x 1]           upper bounds on rho
%
%   Output:
%     fd struct with fields:
%       elemIds          [nTest x 1]  tested rho indices
%       rows             struct array  per-element results
%       maxRelativeError scalar
%       meanRelativeError scalar

    if nargin < 13
        useParallel = false;
    end

    H = numel(rho);
    % Avoid subtracting nearly identical large-solve objectives for weak sensitivities.
    targetObjectiveChange = 1.0e-8;
    candidateIds = find(xmax(:) > xmin(:) + eps);
    nTest = min(nTest, numel(candidateIds));
    elemIds = candidateIds(randperm(numel(candidateIds), nTest));
    rows = repmat(struct('elemId', 0, 'analytic', 0, 'finiteDifference', 0, ...
        'absoluteError', 0, 'relativeError', 0, 'step', 0), nTest, 1);
    relErrors = zeros(nTest, 1);

    for i = 1:nTest
        e = elemIds(i);
        maxCentralStep = min(0.49 * (xmax(e) - rho(e)), 0.49 * (rho(e) - xmin(e)));
        hUse = min(maxCentralStep, max(h, targetObjectiveChange / max(abs(dJdrho(e)), eps)));
        assert(hUse > 0, 'Finite-difference rho(%d) has no free perturbation range.', e);

        rhop = rho;
        rhom = rho;
        rhop(e) = rhop(e) + hUse;
        rhom(e) = rhom(e) - hUse;

        xp = model.segmentToArm(rhop);
        xm = model.segmentToArm(rhom);

        [Jp, ~] = evaluateObjectiveOnly(analyses, xp, penal, pAgg, weights, C0, useParallel);
        [Jm, ~] = evaluateObjectiveOnly(analyses, xm, penal, pAgg, weights, C0, useParallel);
        fdGrad = (Jp - Jm) / (rhop(e) - rhom(e));
        absErr = abs(fdGrad - dJdrho(e));
        relErr = absErr / max([abs(fdGrad), abs(dJdrho(e)), eps]);

        rows(i).elemId = e;
        rows(i).analytic = dJdrho(e);
        rows(i).finiteDifference = fdGrad;
        rows(i).absoluteError = absErr;
        rows(i).relativeError = relErr;
        rows(i).step = hUse;
        relErrors(i) = relErr;

        fprintf('  rho(%7d): analytic=% .6e  FD=% .6e  relErr=%.3e  h=%.1e\n', ...
            e, dJdrho(e), fdGrad, relErr, hUse);
    end

    fd.elemIds = elemIds;
    fd.rows = rows;
    fd.maxRelativeError = max(relErrors);
    fd.meanRelativeError = mean(relErrors);
end
