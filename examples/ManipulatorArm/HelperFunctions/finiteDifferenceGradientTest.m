function fd = finiteDifferenceGradientTest(analyses, x, penal, pAgg, weights, C0, grad, nTest, h, xmin, xmax, useParallel)
    if nargin < 12
        useParallel = false;
    end
    n = numel(x);
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
        maxCentralStep = min(0.49 * (xmax(e) - x(e)), 0.49 * (x(e) - xmin(e)));
        hUse = min(maxCentralStep, max(h, targetObjectiveChange / max(abs(grad(e)), eps)));
        assert(hUse > 0, 'Finite-difference element %d has no free perturbation range.', e);

        xp = x;
        xm = x;
        xp(e) = xp(e) + hUse;
        xm(e) = xm(e) - hUse;

        [Jp, ~] = evaluateObjectiveOnly(analyses, xp, penal, pAgg, weights, C0, useParallel);
        [Jm, ~] = evaluateObjectiveOnly(analyses, xm, penal, pAgg, weights, C0, useParallel);
        fdGrad = (Jp - Jm) / (xp(e) - xm(e));
        absErr = abs(fdGrad - grad(e));
        relErr = absErr / max([abs(fdGrad), abs(grad(e)), eps]);

        rows(i).elemId = e;
        rows(i).analytic = grad(e);
        rows(i).finiteDifference = fdGrad;
        rows(i).absoluteError = absErr;
        rows(i).relativeError = relErr;
        rows(i).step = hUse;
        relErrors(i) = relErr;

        fprintf('  elem %7d: analytic=% .6e  FD=% .6e  relErr=%.3e  h=%.1e\n', ...
            e, grad(e), fdGrad, relErr, hUse);
    end

    fd.elemIds = elemIds;
    fd.rows = rows;
    fd.maxRelativeError = max(relErrors);
    fd.meanRelativeError = mean(relErrors);
end
