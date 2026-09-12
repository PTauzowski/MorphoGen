function configs = extractExtremalConfigs(vN, vTy, vTz, vMs, vMy, vMz, samples)
% EXTRACTEXTREMALCONFIGS  Find statistically extremal arm load configurations.
%
%   configs = extractExtremalConfigs(vN, vTy, vTz, vMs, vMy, vMz, samples)
%
%   Searches a grid of joint-angle samples for the configuration that
%   produces the worst-case internal-force envelope in each of four
%   structural quantities.  The worst-case scalar for sample k is:
%
%       bending:  max over element ends and segments of hypot(My, Mz)
%       torsion:  max over element ends and segments of |Ms|
%       shear:    max over element ends and segments of hypot(Ty, Tz)
%       tension:  max over element ends and segments of N (signed, tensile)
%
%   Under the arm symmetry  forces(-beta) = -forces(beta),  the
%   configuration that minimises a quantity is the exact mirror of the
%   maximiser.  The returned pairs enforce this:
%       min_betas = -max_betas
%
%   Inputs:
%     vN, vTy, vTz, vMs, vMy, vMz  [nSamples x 2 x nElems] internal forces
%     samples                        [nSamples x 7]  joint angles (degrees)
%
%   Output: configs  {7 x 1} cell array of structs with fields
%     .name   string  configuration name
%     .label  string  display label
%     .betas  [1 x 7] joint angles in degrees
%
%   Row order: min_bending, min_torsion, min_shear,
%              max_bending, max_torsion, max_shear, max_tension.

    % Bending resultant envelope: max over both element ends (dim2) and all
    % segments (dim3), then pick sample index.
    Mb_env = max(max(hypot(vMy, vMz), [], 2), [], 3);
    [~, iBend]  = max(Mb_env);
    bBend  = samples(iBend,  :);

    % Torsion envelope
    Ms_env = max(max(abs(vMs), [], 2), [], 3);
    [~, iTors]  = max(Ms_env);
    bTors  = samples(iTors,  :);

    % Shear resultant envelope
    T_env  = max(max(hypot(vTy, vTz), [], 2), [], 3);
    [~, iShear] = max(T_env);
    bShear = samples(iShear, :);

    % Max axial tension (signed N, positive = tensile)
    N_env  = max(max(vN, [], 2), [], 3);
    [~, iTens]  = max(N_env);
    bTens  = samples(iTens,  :);

    configs = {
        struct('name', 'min_bending', 'label', 'Min M_b', 'betas', -bBend);
        struct('name', 'min_torsion', 'label', 'Min M_s', 'betas', -bTors);
        struct('name', 'min_shear',   'label', 'Min T',   'betas', -bShear);
        struct('name', 'max_bending', 'label', 'Max M_b', 'betas',  bBend);
        struct('name', 'max_torsion', 'label', 'Max M_s', 'betas',  bTors);
        struct('name', 'max_shear',   'label', 'Max T',   'betas',  bShear);
        struct('name', 'max_tension', 'label', 'Max N',   'betas',  bTens);
    };

    fprintf('[extractExtremalConfigs] max_bending  betas = %s\n', mat2str(bBend,  4));
    fprintf('[extractExtremalConfigs] max_torsion  betas = %s\n', mat2str(bTors,  4));
    fprintf('[extractExtremalConfigs] max_shear    betas = %s\n', mat2str(bShear, 4));
    fprintf('[extractExtremalConfigs] max_tension  betas = %s\n', mat2str(bTens,  4));
end
