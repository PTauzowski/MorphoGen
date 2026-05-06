function [originRef, originNames, confidence] = classifyCurveFamilyOrigin(fieldsRef, mixedRatio)
% classifyCurveFamilyOrigin  Assign each reference element to dominant family.
%
% originRef = 0 marks mixed/ambiguous elements. Positive labels correspond to
% originNames entries.

    if nargin < 2
        mixedRatio = 0.85;
    end

    originNames = ["helixPlus"; "helixMinus"; "axial"; "bending"; "ring"; "jointRing"];
    n = numel(fieldsRef.envelope);
    values = zeros(n, numel(originNames));
    for i = 1:numel(originNames)
        name = char(originNames(i));
        if isfield(fieldsRef, name)
            values(:, i) = fieldsRef.(name)(:);
        end
    end

    [sortedValues, sortedIds] = sort(values, 2, 'descend');
    best = sortedValues(:, 1);
    second = sortedValues(:, 2);

    originRef = sortedIds(:, 1);
    ambiguous = best <= eps | (second ./ max(best, eps)) > mixedRatio;
    originRef(ambiguous) = 0;
    confidence = 1 - second ./ max(best, eps);
    confidence(best <= eps) = 0;
end
