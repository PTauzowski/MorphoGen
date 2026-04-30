function fd = finiteDifferenceGradientTestLinked(model, analyses, rho, penal, pAgg, weights, C0, dJdrho, nTest, h, xmin, xmax)
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

    H = numel(rho);
    nTest = min(nTest, H);
    elemIds = randperm(H, nTest)';
    rows = repmat(struct('elemId', 0, 'analytic', 0, 'finiteDifference', 0, ...
        'relativeError', 0), nTest, 1);
    relErrors = zeros(nTest, 1);

    for i = 1:nTest
        e = elemIds(i);
        hUse = min([h, 0.49 * (xmax(e) - rho(e)), 0.49 * (rho(e) - xmin(e))]);
        if hUse <= 0
            hUse = h;
        end

        rhop = rho;
        rhom = rho;
        rhop(e) = min(xmax(e), rhop(e) + hUse);
        rhom(e) = max(xmin(e), rhom(e) - hUse);

        xp = model.segmentToArm(rhop);
        xm = model.segmentToArm(rhom);

        [Jp, ~] = evaluateObjectiveOnly(analyses, xp, penal, pAgg, weights, C0);
        [Jm, ~] = evaluateObjectiveOnly(analyses, xm, penal, pAgg, weights, C0);
        fdGrad = (Jp - Jm) / (rhop(e) - rhom(e));
        relErr = abs(fdGrad - dJdrho(e)) / max([abs(fdGrad), abs(dJdrho(e)), eps]);

        rows(i).elemId = e;
        rows(i).analytic = dJdrho(e);
        rows(i).finiteDifference = fdGrad;
        rows(i).relativeError = relErr;
        relErrors(i) = relErr;

        fprintf('  rho(%7d): analytic=% .6e  FD=% .6e  relErr=%.3e\n', ...
            e, dJdrho(e), fdGrad, relErr);
    end

    fd.elemIds = elemIds;
    fd.rows = rows;
    fd.maxRelativeError = max(relErrors);
    fd.meanRelativeError = mean(relErrors);
end
