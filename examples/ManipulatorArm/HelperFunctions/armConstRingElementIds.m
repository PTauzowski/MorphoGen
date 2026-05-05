function const_elems = armConstRingElementIds(model, arm, designSpace)
% armConstRingElementIds  Constant ring element IDs for arm topologies.
%
%   const_elems = armConstRingElementIds(model, arm, "reference") returns
%   IDs for a ReferenceModuleSolidModel.
%
%   const_elems = armConstRingElementIds(model, arm, "linked") returns IDs
%   in reference half-segment design space.
%
%   const_elems = armConstRingElementIds(model, arm, "full") returns IDs in
%   full-arm element space.

    if nargin < 3
        designSpace = "full";
    end
    designSpace = string(designSpace);

    useEndRing = isfield(arm, 'constEndRing') && arm.constEndRing;
    useMiddleRing = isfield(arm, 'constMiddleRing') && arm.constMiddleRing;

    if ~(useEndRing || useMiddleRing)
        const_elems = zeros(0, 1);
        return;
    end

    if isprop(model, 'loaded_node_ids') && isprop(model, 'fixedFaceSelector')
        const_elems = referenceModuleConstElems(model, useEndRing, useMiddleRing);
        return;
    end

    H = model.halfSegmentNelems;
    sliceElems = localSliceElementCount(arm);
    assert(sliceElems <= H, ...
        'Ring slice size %d exceeds half-segment element count %d.', sliceElems, H);

    halfIds = halfSegmentRingIds(H, sliceElems, useEndRing, useMiddleRing);
    halfIds = unique(halfIds(:));

    switch designSpace
        case "linked"
            const_elems = halfIds;
        case "full"
            nElems = model.analysis.getTotalElemsNumber();
            nCopies = nElems / H;
            assert(abs(nCopies - round(nCopies)) < eps, ...
                'Full-arm element count %d is not an integer multiple of H=%d.', nElems, H);
            const_elems = fullArmPhysicalRingIds(H, round(nCopies), sliceElems, ...
                useEndRing, useMiddleRing);
        otherwise
            error('Unknown designSpace "%s". Use "reference", "linked", or "full".', designSpace);
    end
end

function const_elems = referenceModuleConstElems(model, useEndRing, useMiddleRing)
    nodeIds = zeros(0, 1);
    if useEndRing
        nodeIds = [nodeIds; model.loaded_node_ids(:)]; %#ok<AGROW>
    end
    if useMiddleRing
        fixedNodeIds = find(model.fixedFaceSelector.select(model.mesh.nodes));
        nodeIds = [nodeIds; fixedNodeIds(:)]; %#ok<AGROW>
    end
    const_elems = find(any(ismember(model.mesh.elems, unique(nodeIds)), 2));
    const_elems = unique(const_elems(:));
end

function sliceElems = localSliceElementCount(arm)
    resTh = max(1, round(arm.res_th));
    wallThickness = arm.R - arm.r;
    resCirc = max(3, round(2 * pi * arm.R / wallThickness * resTh));
    sliceElems = resTh * resCirc;
end

function ids = halfSegmentRingIds(H, sliceElems, useEndRing, useMiddleRing)
    ids = zeros(0, 1);
    if useEndRing
        ids = [ids; (1:sliceElems)']; %#ok<AGROW>
    end
    if useMiddleRing
        ids = [ids; (H - sliceElems + 1:H)']; %#ok<AGROW>
    end
end

function ids = fullArmPhysicalRingIds(H, nCopies, sliceElems, useEndRing, useMiddleRing)
    ids = zeros(0, 1);

    % Boundary numbering for 12 half-segments:
    %   b=0,2,4,... are segment ends; b=1,3,5,... are segment middles.
    % Use one axial element layer per physical boundary, not both adjacent
    % layers around an internal boundary.
    if useEndRing
        ids = [ids; (1:sliceElems)']; %#ok<AGROW> % root/end boundary b=0
        for b = 2:2:nCopies
            blockStart = (b - 1) * H + 1;
            blockEnd = blockStart + H - 1;
            ids = [ids; (blockEnd - sliceElems + 1:blockEnd)']; %#ok<AGROW>
        end
    end

    if useMiddleRing
        for b = 1:2:nCopies-1
            blockStart = (b - 1) * H + 1;
            blockEnd = blockStart + H - 1;
            ids = [ids; (blockEnd - sliceElems + 1:blockEnd)']; %#ok<AGROW>
        end
    end

    ids = unique(ids(:));
end
