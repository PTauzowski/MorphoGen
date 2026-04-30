function [C, dC] = computeComplianceAndGradient(analysis, x, penal)
    nElemsTotal = analysis.getTotalElemsNumber();
    assert(numel(x) == nElemsTotal, ...
        'Density vector length %d does not match analysis element count %d.', ...
        numel(x), nElemsTotal);

    x = x(:);
    qfem = analysis.solveWeighted(x .^ penal);
    q = analysis.fromFEMVector(qfem(:, 1));
    P = analysis.Pfem(:, 1);
    C = P' * qfem(:, 1);

    if analysis.isConst
        stiffnessFunction = 'computeStifnessMatrixConst';
    else
        stiffnessFunction = 'computeStifnessMatrix';
    end

    dC = zeros(nElemsTotal, 1);
    elemOffset = 0;
    xOnes = ones(nElemsTotal, 1);
    elemIndices = analysis.getElemIndices();

    for i = 1:numel(analysis.felems)
        fe = analysis.felems{i};
        elemIds = elemIndices{i};
        nelems = size(fe.elems, 1);
        nnodes = size(fe.elems, 2);
        ndofs = size(fe.ndofs, 2);
        dim = nnodes * ndofs;
        K0 = reshape(fe.(stiffnessFunction)(analysis.mesh.nodes, xOnes(elemIds)), ...
            dim, dim, nelems);
        qelems = fe.createElemSolutionVectors(q);

        for e = 1:nelems
            globalElem = elemOffset + e;
            elemEnergy = qelems(:, e)' * K0(:, :, e) * qelems(:, e);
            dC(globalElem) = -penal * x(globalElem)^(penal - 1) * elemEnergy;
        end
        elemOffset = elemOffset + nelems;
    end
end
