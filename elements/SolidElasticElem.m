function K = computeStifnessMatrix(obj, nodes, el_idx)

    % ---- Optional argument handling ----
    if nargin < 3 || isempty(el_idx)
        el_idx = [];
    end

    nelems = size(obj.elems,1);
    if ~isempty(el_idx)
        nelems = numel(el_idx);
    end

    nnodes = size(obj.elems,2);
    ndofs  = size(obj.ndofs,2);
    dim    = nnodes * ndofs;

    integrator = obj.sf.createIntegrator();
    nip = size(integrator.points,1);

    dN = permute( ...
        repmat(obj.sf.computeGradient(integrator.points),[1,1,1,nelems]), ...
        [2,1,4,3]);

    [~, J1, detJ] = obj.computeJacobian(nodes, dN, el_idx);
    dNx = pagemtimes(J1, dN);

    B = obj.computeStrainDerivativesMatrix(dNx, nip, el_idx);

    h = repelem(obj.props.h, nelems, 1);

    K = obj.computeElementMatrices(h, integrator.weights, detJ, B, obj.mat.D);
end
