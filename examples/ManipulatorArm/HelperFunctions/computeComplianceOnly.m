function C = computeComplianceOnly(analysis, x, penal)
    nElemsTotal = analysis.getTotalElemsNumber();
    assert(numel(x) == nElemsTotal, ...
        'Density vector length %d does not match analysis element count %d.', ...
        numel(x), nElemsTotal);

    qfem = analysis.solveWeighted(x(:) .^ penal);
    P = analysis.Pfem(:, 1);
    C = P' * qfem(:, 1);
end
