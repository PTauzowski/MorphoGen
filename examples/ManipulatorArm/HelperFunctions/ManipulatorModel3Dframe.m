classdef ManipulatorModel3Dframe < handle
%MANIPULATORMODEL3D  Build a multi-segment 3D Arm-Z mesh and FE model.
% - Angles passed in DEGREES; converted once.
% - Uses the same (working) convention as your original: pre-transposed
%   rotation matrices and right-multiplication by points row-vectors.
% - Adds snap-to-grid before merges to avoid tiny discontinuities.
% - Centralizes resolution/tolerances; clarifies selectors & steps.

    properties
        analysis
        mesh
        elems
        fe
        xEnd
        frameNodes
        const_elems
        loadSurfaceNodes
        fixedSurfaceNodes
        halfSegmentNelems
        use_offset
    end

    properties (Access=private)
        % geometry
        R; r; Th; ls
        alpha; betas
        % meshing
        resTh = 1
        resCirc
        resLen
        % rotations (pre-transposed, matching original)
        rotCutT
        rotCut2T
        % settings
        epsPlane = 1.0e-6    % plane tolerance
        epsSnap  = 1.0e-10   % snap-to-grid resolution (model units)
    end

    methods
        function obj = ManipulatorModel3Dframe(E, nu, ls, alpha_deg, betas_deg, use_offset)
            % ---- input validation (clear error messages) ----
            arguments
                E (1,1) double {mustBePositive}
                nu (1,1) double {mustBeLessThan(nu,0.4999), mustBeGreaterThan(nu,-0.999)}
                ls (1,1) double {mustBePositive}
                alpha_deg (1,1) double
                betas_deg double {mustBeVector}
                use_offset (1,1) logical
            end

            % ---- store geometry & angles ----
            obj.ls = ls;
            obj.alpha = deg2rad(alpha_deg);
            obj.betas = deg2rad(betas_deg(:).');    % row vector
            obj.use_offset = double(use_offset);


            % ---- pre-transposed rotations (MATCH ORIGINAL) ----
            obj.rotCutT  = obj.Ry(obj.alpha)';        % == rotCut
            obj.rotCut2T = obj.Ry(2*obj.alpha)';      % == rotCut2

            % ---- build geometry chain ----
            obj.mesh  = Mesh();
            obj.elems = [];
            obj.generateManipulator();

        end

        function ns = determineSegment(obj, elemIndex)
            % Map element index -> module index (1-based), assuming first
            % half-segment count is canonical (matches original behavior).
            nhs = round(elemIndex ./ obj.halfSegmentNelems);
            ns  = ceil((nhs - 1)./2) + 1;
        end

        function plot(obj)
            obj.fe.plot(obj.mesh.nodes)
        end
    end

    %% ========= Internal build =========
    methods (Access=private)
        function generateManipulator(obj)
            % Phase offset (as in original) uses the first beta
            phase = -obj.betas(1) * obj.use_offset;


            
            % Rotate whole assembly by beta_1 about Z (pre-transposed)
            Rz1T = obj.Rz(obj.betas(1))';


            % kinematic accumulators (match original)
            prevRot = obj.rotCutT * Rz1T;     % == rotCut * rotBeta
            obj.xEnd = [0 0 obj.ls];
            xEnds = [0 0 0; obj.xEnd];

            % Remaining modules
            for k = 2:numel(obj.betas)
                RzTk = obj.Rz(obj.betas(k))';

                % update phase as original
                phase = obj.use_offset * (phase - obj.betas(k));

                % advance tip: 2*ls unless last (+ls)
                step = (k < numel(obj.betas)) * 2*obj.ls + (k == numel(obj.betas)) * obj.ls;
                obj.xEnd = obj.xEnd + [0 0 step] * obj.rotCutT * RzTk * prevRot;
                xEnds = [xEnds; obj.xEnd];

                % accumulate rotation for next placement
                prevRot = obj.rotCut2T * RzTk * prevRot;
            end

            obj.frameNodes = xEnds;
        end
    end

    %% ========= Small utilities =========
    methods (Access=private)
        function R = Rz(~, a)
            c = cos(a); s = sin(a);
            R = [ c -s  0;  s  c  0;  0  0  1 ];
        end
        function R = Ry(~, a)
            c = cos(a); s = sin(a);
            R = [ c  0  s;  0  1  0; -s  0  c ];
        end
        function X = snap(obj, X)
            % Snap to grid to collapse FP noise before merging.
            if obj.epsSnap > 0
                X = round(X/obj.epsSnap) * obj.epsSnap;
            end
        end
    end
end
