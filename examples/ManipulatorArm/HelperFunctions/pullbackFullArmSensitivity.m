function dJdrho = pullbackFullArmSensitivity(dJdx, H, nElems)
% PULLBACKFULLARMSENSITIVITY  Pull back full-arm sensitivity to rho-space
%   via the exact chain rule of the segmentToArm expansion.
%
%   dJdrho = pullbackFullArmSensitivity(dJdx, H, nElems)
%
%   Inputs:
%     dJdx   [nElems x 1]  sensitivity of objective w.r.t. every full-arm element
%     H      scalar        number of elements in one reference half-segment
%     nElems scalar        full-arm element count  (= 2 * nArms * H)
%
%   Output:
%     dJdrho [H x 1]  sensitivity w.r.t. the reference half-segment densities
%
%   Derivation:
%     x_arm = segmentToArm(rho) expands rho into nCopies = 2*nArms contiguous
%     blocks of H elements each.  Odd-numbered blocks are normal copies
%     (x_arm((c-1)*H+p) = rho(p)); even-numbered are flipped
%     (x_arm((c-1)*H+p) = rho(H+1-p)).
%
%     By the chain rule:
%       dJ/drho(e) = sum_{e_global : x_arm(e_global) = rho(e)} dJ/dx(e_global)
%
%     For normal copy c at local position p=e:   contributes dJdx((c-1)*H+e)
%     For flipped copy c at local position p=H+1-e: contributes dJdx((c-1)*H+(H+1-e))
%
%     Vectorized: reshape dJdx into [H x nCopies], sum normal columns
%     directly and sum flipped columns after flipping along the row dimension.

    assert(mod(nElems, H) == 0, 'nElems (%d) must be divisible by H (%d).', nElems, H);
    nCopies = nElems / H;

    % Reshape: column c holds sensitivities for copy c, rows are local positions
    S = reshape(dJdx(:), H, nCopies);   % [H x nCopies]

    normalCols  = S(:, 1:2:end);        % [H x nArms]  odd copies
    flippedCols = S(:, 2:2:end);        % [H x nArms]  even copies

    % Normal copy contribution: local position p → rho(p), add directly
    % Flipped copy contribution: local position p → rho(H+1-p), flip before sum
    dJdrho = sum(normalCols, 2) + sum(flip(flippedCols, 1), 2);
end
