function ais_rho = pullbackFullArmAverageIntensity(ais_full, H, nElems)
% PULLBACKFULLARMAVERAGEINTENSITY  Pull back full-arm element intensities
%   to rho-space by averaging over all copies that share the same rho entry.
%
%   ais_rho = pullbackFullArmAverageIntensity(ais_full, H, nElems)
%
%   Inputs:
%     ais_full  [nElems x 1]  per-element intensity on the full arm
%                              (already max-over-configs, filtered)
%     H         scalar        half-segment element count (reference module)
%     nElems    scalar        full-arm element count (= 2*nArms*H)
%
%   Output:
%     ais_rho  [H x 1]  average intensity of rho(e_ref) across all its copies
%
%   The segmentToArm expansion produces nCopies = nElems/H contiguous blocks:
%     odd blocks  → normal copy:  x_arm((c-1)*H+p) = rho(p)
%     even blocks → flipped copy: x_arm((c-1)*H+p) = rho(H+1-p)
%
%   For rho(e_ref), the contributing full-arm elements are:
%     - normal copies:  global index (c-1)*H + e_ref
%     - flipped copies: global index (c-1)*H + (H+1-e_ref)
%
%   The pullback computes the mean intensity over all nCopies contributing
%   elements so that the removal criterion reflects the average load carried
%   by each reference module element across all arm positions.

    assert(mod(nElems, H) == 0, 'nElems (%d) must be divisible by H (%d).', nElems, H);
    nCopies = nElems / H;

    % Reshape to [H x nCopies]: column c = intensities for copy c
    S = reshape(ais_full(:), H, nCopies);

    normalCols  = S(:, 1:2:end);   % [H x nArms] — odd copies (normal)
    flippedCols = S(:, 2:2:end);   % [H x nArms] — even copies (flipped)

    % Normal copy at position p contributes to rho(p)  → add directly
    % Flipped copy at position p contributes to rho(H+1-p) → flip before averaging
    ais_rho = mean([normalCols, flip(flippedCols, 1)], 2);
end
