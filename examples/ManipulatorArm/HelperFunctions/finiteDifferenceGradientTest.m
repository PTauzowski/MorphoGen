function fd = finiteDifferenceGradientTest(analyses, x, penal, pAgg, weights, C0, grad, nTest, h, xmin, xmax)
    n = numel(x);
    nTest = min(nTest, n);
    elemIds = randperm(n, nTest)';
    rows = repmat(struct('elemId', 0, 'analytic', 0, 'finiteDifference', 0, ...
        'relativeError', 0), nTest, 1);
    relErrors = zeros(nTest, 1);

    for i = 1:nTest
        e = elemIds(i);
        hUse = min([h, 0.49 * (xmax(e) - x(e)), 0.49 * (x(e) - xmin(e))]);
        if hUse <= 0
            hUse = h;
        end

        xp = x;
        xm = x;
        xp(e) = min(xmax(e), xp(e) + hUse);
        xm(e) = max(xmin(e), xm(e) - hUse);

        [Jp, ~] = evaluateObjectiveOnly(analyses, xp, penal, pAgg, weights, C0);
        [Jm, ~] = evaluateObjectiveOnly(analyses, xm, penal, pAgg, weights, C0);
        fdGrad = (Jp - Jm) / (xp(e) - xm(e));
        relErr = abs(fdGrad - grad(e)) / max([abs(fdGrad), abs(grad(e)), eps]);

        rows(i).elemId = e;
        rows(i).analytic = grad(e);
        rows(i).finiteDifference = fdGrad;
        rows(i).relativeError = relErr;
        relErrors(i) = relErr;

        fprintf('  elem %7d: analytic=% .6e  FD=% .6e  relErr=%.3e\n', ...
            e, grad(e), fdGrad, relErr);
    end

    fd.elemIds = elemIds;
    fd.rows = rows;
    fd.maxRelativeError = max(relErrors);
    fd.meanRelativeError = mean(relErrors);
end
