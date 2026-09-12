function labelsFull = expandReferenceLabelsToArm(model, labelsRef)
% expandReferenceLabelsToArm  Expand reference half-segment labels like segmentToArm.

    labelsRef = labelsRef(:);
    H = numel(labelsRef);
    nElems = model.analysis.getTotalElemsNumber();
    nCopies = nElems / H;
    assert(abs(nCopies - round(nCopies)) < eps, ...
        'Full-arm element count is not an integer multiple of reference labels.');
    nCopies = round(nCopies);

    labelsFull = labelsRef;
    for k = 1:(nCopies/2 - 1)
        labelsFull = [labelsFull; flip(labelsRef); labelsRef]; %#ok<AGROW>
    end
    labelsFull = [labelsFull; flip(labelsRef)];
end
