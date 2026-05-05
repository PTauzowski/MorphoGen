function assertLinkedArmLayoutCompatible(model, referenceModel, configName)
% assertLinkedArmLayoutCompatible  Validate linked design-space compatibility.
%
% Rotation-aware linked models (use_offset=false) intentionally generate each
% half-segment in its local frame. A nonzero joint rotation can shift which
% circumferential node IDs are shared at an interface, so raw mesh.elems need
% not be identical across configurations. The SIMP linkage only requires the
% same ordered half-segment element blocks and conforming interfaces.

    if nargin < 3
        configName = "configuration";
    end
    configName = char(string(configName));

    H = model.halfSegmentNelems;
    Href = referenceModel.halfSegmentNelems;
    assert(H == Href, ...
        'Half-segment element count mismatch in %s: got %d, expected %d.', ...
        configName, H, Href);

    nElems = model.analysis.getTotalElemsNumber();
    nElemsRef = referenceModel.analysis.getTotalElemsNumber();
    assert(nElems == nElemsRef, ...
        'Element-count mismatch in %s: got %d, expected %d.', ...
        configName, nElems, nElemsRef);

    nCopies = nElems / H;
    nCopiesRef = nElemsRef / Href;
    assert(abs(nCopies - round(nCopies)) < eps && abs(nCopiesRef - round(nCopiesRef)) < eps, ...
        'Full-arm element count is not an integer multiple of H in %s.', configName);
    assert(round(nCopies) == round(nCopiesRef), ...
        'Half-segment copy-count mismatch in %s: got %d, expected %d.', ...
        configName, round(nCopies), round(nCopiesRef));

    assert(model.resTh == referenceModel.resTh && ...
           model.resCirc == referenceModel.resCirc && ...
           model.resLen == referenceModel.resLen, ...
        'Mesh resolution mismatch in %s.', configName);

    expectedShared = 2 * model.resCirc * model.resTh;
    for c = 1:round(nCopies)-1
        elemsLeft = (c - 1) * H + (1:H);
        elemsRight = c * H + (1:H);
        shared = intersect(unique(model.mesh.elems(elemsLeft, :)), ...
                           unique(model.mesh.elems(elemsRight, :)));
        assert(numel(shared) == expectedShared, ...
            'Interface %d-%d in %s has %d shared nodes, expected %d.', ...
            c, c + 1, configName, numel(shared), expectedShared);
    end
end
