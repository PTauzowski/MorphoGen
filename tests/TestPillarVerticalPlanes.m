classdef TestPillarVerticalPlanes < matlab.unittest.TestCase
    % TestPillarVerticalPlanes  Verifies that PillarModel generates nodes at
    % every canonical z-plane dictated by the interface-layer rule.
    %
    % KEY INVARIANT: for every pair of adjacent nominal layers, an interface
    % slab of thickness int_th must exist.  The canonical z-planes are the
    % boundaries of every slab (nominal or interface).  Every region of the
    % mesh that spans a canonical z-plane must contain nodes exactly at that
    % z-level.
    %
    % This test is designed to FAIL on the pre-fix code (which was missing
    % z = +int_th/2 and z = z_pil1_top+int_th in the bank ring3D/transition
    % blocks) and to PASS after the surgical fix.

    methods (TestClassSetup)
        function addPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            dirs = {'math', 'mesh', 'analysis', 'elements', ...
                    fullfile('examples', 'Pillar')};
            for d = dirs
                testCase.applyFixture( ...
                    matlab.unittest.fixtures.PathFixture(fullfile(root, d{1})));
            end
        end
    end

    methods (Test)
        function testBankAndCoreCanonicalZPlanes(testCase)
            % --- small model for fast execution ---
            top_R         = 5;
            pillar_layers = [4, 8];
            pillar_res    = [1, 2];
            ground_layers = [10, 6];
            ground_res    = [3, 2];
            pillar_chem   = [0.08, 0.08];
            ground_chem   = [0, 1];
            int_th        = 1.0;
            z_offset      = 0;

            opts.nrXY      = 2;    % reduced resolution: fast test
            opts.bank_hres = 2;
            opts.nThetaSeg = 4;

            sf    = ShapeFunctionH27();
            model = PillarModel(top_R, pillar_layers, ground_layers, ...
                                pillar_res, ground_res, ...
                                pillar_chem, ground_chem, ...
                                int_th, sf, z_offset, opts);

            nodes = model.mesh.nodes;
            tolZ  = 1e-5;

            % ----------------------------------------------------------------
            % 1.  Compute canonical z-planes analytically (same arithmetic as
            %     addLayeredQuarterCylinder / addLayeredPipe3D) and verify each
            %     plane has at least one mesh node.  This catches any globally
            %     missing plane.
            % ----------------------------------------------------------------
            zCanon = computeCanonicalZ(int_th, ground_layers, pillar_layers);
            zAll   = nodes(:, 3);

            for j = 1:numel(zCanon)
                testCase.verifyTrue( ...
                    any(abs(zAll - zCanon(j)) <= tolZ), ...
                    sprintf('No mesh node found at canonical z = %.6f', zCanon(j)));
            end

            % ----------------------------------------------------------------
            % 2.  Check the bank outer wall (r ~ Rbank).
            %
            %     The outer radius of the bank ring3D blocks is Rbank (constant
            %     for all z above -int_th/2).  The addLayeredPipe3D for the
            %     ground also uses Rbank as the outer radius.  After weldNodes
            %     these share nodes, so this column spans the full model height.
            %
            %     The THREE planes listed in `critical` were absent in the
            %     buggy code:
            %       z_seam_top   = +int_th/2    (top of global seam interface)
            %       z_pil1_top   = pillar_layers(1) - int_th/2  (effective top)
            %       z_pil1_int_t = z_pil1_top + int_th          (internal intf top)
            % ----------------------------------------------------------------
            pillar_height = sum(pillar_layers);
            Rbase = top_R + pillar_height * tan(5 * pi / 180);
            Rout  = 1.2 * Rbase;
            Rbank = 1.5 * Rout;

            r    = hypot(nodes(:, 1), nodes(:, 2));
            tolR = max(1e-4, Rbank * 5e-4);
            bankZ = nodes(abs(r - Rbank) <= tolR, 3);

            testCase.assertFalse(isempty(bankZ), ...
                'No nodes found near r = Rbank; check geometry parameters.');

            z_seam_top   = int_th / 2;
            z_pil1_top   = pillar_layers(1) - int_th / 2;
            z_pil1_int_t = z_pil1_top + int_th;
            critical = [z_seam_top; z_pil1_top; z_pil1_int_t];

            labels = {'z_seam_top (+int_th/2)', ...
                      'z_pil1_top (pil_layers(1)-int_th/2)', ...
                      'z_pil1_int_top (z_pil1_top+int_th)'};
            for j = 1:numel(critical)
                testCase.verifyTrue( ...
                    any(abs(bankZ - critical(j)) <= tolZ), ...
                    sprintf('Bank (r~Rbank): missing %s = %.6f', ...
                            labels{j}, critical(j)));
            end

            % ----------------------------------------------------------------
            % 3.  Check the transition outer boundary (x ~ Rtile or y ~ Rtile).
            %
            %     addLayeredTransitionCircleToSquare builds a square-tile region
            %     whose outer edge is a square of half-width Rtile.  The same
            %     critical planes must appear there.
            % ----------------------------------------------------------------
            Rtile = 1.2 * Rbank;
            tolT  = max(1e-4, Rtile * 5e-4);

            onOuter = (abs(nodes(:,1) - Rtile) <= tolT) | ...
                      (abs(nodes(:,2) - Rtile) <= tolT);
            % Restrict to the upper-bank z range (above the ground top)
            aboveGround = nodes(:, 3) >= -int_th / 2 - tolZ;
            transZ = nodes(onOuter & aboveGround, 3);

            % The outer tile boundary is only created if nThetaSeg > 0; skip
            % gracefully if no nodes are found (degenerate resolution).
            if ~isempty(transZ)
                for j = 1:numel(critical)
                    testCase.verifyTrue( ...
                        any(abs(transZ - critical(j)) <= tolZ), ...
                        sprintf('Transition (x/y~Rtile): missing %s = %.6f', ...
                                labels{j}, critical(j)));
                end
            end
        end
    end

end

% =========================================================================
%  Local helper functions
% =========================================================================

function zPlanes = computeCanonicalZ(int_th, ground_layers, pillar_layers)
    % Canonical element-face z-planes for a PillarModel with z_offset = 0.
    % Mirrors the exact h_eff arithmetic used by addLayeredQuarterCylinder
    % and addLayeredPipe3D — NOT derived from unique(mesh.nodes(:,3)).
    %
    % Ground stack:  starts at z = -sum(ground_layers),
    %                last nominal layer trimmed by int_th/2 before entering
    %                the global seam, ends at z = -int_th/2.
    % Global seam:   z in [-int_th/2, +int_th/2].
    % Pillar stack:  starts at z = +int_th/2,
    %                first nominal layer trimmed by int_th/2,
    %                ends at z = sum(pillar_layers).

    Lg = numel(ground_layers);
    Lp = numel(pillar_layers);

    % --- ground effective inputs (last layer trimmed for global seam) ---
    ge    = ground_layers(:);
    ge(Lg) = ge(Lg) - int_th / 2;

    % --- internal h_eff for ground (addLayeredPipe3D / addLayeredQuarterCylinder) ---
    hg = ge;
    for k = 1:Lg
        if k > 1,  hg(k) = hg(k) - int_th / 2; end
        if k < Lg, hg(k) = hg(k) - int_th / 2; end
    end

    % --- pillar effective inputs (first layer trimmed for global seam) ---
    pe    = pillar_layers(:);
    pe(1) = pe(1) - int_th / 2;

    % --- internal h_eff for pillar ---
    hp = pe;
    for k = 1:Lp
        if k > 1,  hp(k) = hp(k) - int_th / 2; end
        if k < Lp, hp(k) = hp(k) - int_th / 2; end
    end

    % --- accumulate z-planes ---
    z_bot  = -sum(ground_layers);
    zPlanes = z_bot;
    z = z_bot;

    % Ground layers and their inter-layer interfaces
    for k = 1:Lg
        z = z + hg(k);
        zPlanes(end+1) = z; %#ok<AGROW>
        if k < Lg
            z = z + int_th;       % internal ground interface
            zPlanes(end+1) = z;   %#ok<AGROW>
        end
    end
    % z == -int_th/2 here

    % Global seam top
    z = z + int_th;               % z == +int_th/2
    zPlanes(end+1) = z;           %#ok<AGROW>

    % Pillar layers and their inter-layer interfaces
    for k = 1:Lp
        z = z + hp(k);
        zPlanes(end+1) = z;       %#ok<AGROW>
        if k < Lp
            z = z + int_th;       % internal pillar interface
            zPlanes(end+1) = z;   %#ok<AGROW>
        end
    end
    % z == sum(pillar_layers) here

    zPlanes = sort(unique(zPlanes));
end
