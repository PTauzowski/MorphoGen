classdef LinearElasticityWeighted < LinearElasticity
    % DEPRECATED transitional shim -- D1, removed in Phase 6.
    %
    % D7 folded this class into LinearElasticity. 61 files on this branch still
    % name it and 47 still call solveWeighted, and migrating them is Phase 6's
    % job, so the name survives the merge as a thin subclass that adds nothing.
    % Deleting it here would have broken every one of those call sites inside a
    % merge resolution, which is exactly what Rule 1 forbids bundling.
    methods
        function obj = LinearElasticityWeighted(felems, mesh, isConst)
            obj = obj@LinearElasticity(felems, mesh, isConst);
        end

        function qfem = solveWeighted(obj, x, retainStiffness)
            if nargin < 3, retainStiffness = false; end
            qfem = obj.solve(x, retainStiffness);
        end
    end
end
