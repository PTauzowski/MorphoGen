function map = buildLinkedSensitivityMap(H, nElems)
% BUILDLINKEDSENSITIVITYMAP  Build index map from reference half-segment to
%   full-arm elements under the segmentToArm expansion pattern.
%
%   map = buildLinkedSensitivityMap(H, nElems)
%
%   Inputs:
%     H       scalar  number of elements in one reference half-segment
%     nElems  scalar  total full-arm element count (must be even multiple of H)
%
%   Output:
%     map struct with fields:
%       normalMap   [nArms x H]  map.normalMap(i, e_ref)  = global index of
%                                  the i-th normal copy of rho(e_ref)
%       flippedMap  [nArms x H]  map.flippedMap(i, e_ref) = global index of
%                                  the i-th flipped copy of rho(e_ref)
%       H, nElems, nArms, nCopies  (scalars, for verification)
%
%   The segmentToArm expansion is:
%     [rho; flip(rho); rho; flip(rho); ...]  (2*nArms copies total)
%   Odd-numbered copies are normal; even-numbered copies are flipped.
%
%   For normal copy c (odd) at local position p:  x_arm((c-1)*H+p) = rho(p)
%   For flipped copy c (even) at local position p: x_arm((c-1)*H+p) = rho(H+1-p)
%
%   Therefore perturbing rho(e_ref) changes exactly:
%     - All global elements in map.normalMap(:, e_ref)
%     - All global elements in map.flippedMap(:, e_ref)

    assert(H > 0 && mod(nElems, H) == 0, ...
        'nElems (%d) must be a positive integer multiple of H (%d).', nElems, H);
    nCopies = nElems / H;
    assert(mod(nCopies, 2) == 0, ...
        'nCopies = nElems/H = %d must be even (half-segment pairs).', nCopies);
    nArms = nCopies / 2;

    % Copy-number bases: global start index for copy c is (c-1)*H + 1
    normalCopyNums  = (1:2:nCopies)';    % odd copies  [nArms x 1]
    flippedCopyNums = (2:2:nCopies)';    % even copies [nArms x 1]
    normalBases     = (normalCopyNums  - 1) * H;   % [nArms x 1]
    flippedBases    = (flippedCopyNums - 1) * H;   % [nArms x 1]

    % Normal copy c, position p: global = (c-1)*H + p  → rho(p)
    %   normalMap(i, e_ref) = normalBases(i) + e_ref
    eRef = 1:H;                                            % [1 x H]
    normalMap = bsxfun(@plus, normalBases, eRef);          % [nArms x H]

    % Flipped copy c, position p: global = (c-1)*H + p  → rho(H+1-p)
    %   So rho(e_ref) sits at flipped-copy position p = H+1-e_ref
    %   flippedMap(i, e_ref) = flippedBases(i) + (H+1-e_ref)
    flippedMap = bsxfun(@plus, flippedBases, H + 1 - eRef);  % [nArms x H]

    map.normalMap  = normalMap;
    map.flippedMap = flippedMap;
    map.H          = H;
    map.nElems     = nElems;
    map.nArms      = nArms;
    map.nCopies    = nCopies;
end
