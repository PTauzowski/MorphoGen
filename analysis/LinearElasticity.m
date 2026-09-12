classdef LinearElasticity < FEAnalysis
    % Linear elastic analysis, optionally density-weighted.
    %
    % D7 folds LinearElasticityWeighted into this class. The port runs three
    % ways, because each branch had something the others lacked:
    %
    %   from CAS_Arm      the assembled-stiffness cache and solveAdjointWithLoad,
    %                     which exists to avoid reassembling K for the adjoint
    %                     system, plus saveMatrices
    %   from Vibrations   selfLoadFactor and the self-weight load path
    %   from develop      the prepareRHSVectors() call, which Vibrations had
    %                     commented out -- an omission rather than a decision:
    %                     develop and CAS_Arm both call it and six other
    %                     Vibrations analyses still do. It is the step that
    %                     moves accumulated nodal load into the RHS.

    properties
        isConst;
        selfLoadFactor;
        cachedWeightedX;
        cachedWeightedKvals;
        cachedWeightedFunction;
    end

    methods
        function obj = LinearElasticity(felems, mesh, isConst)
            obj = obj@FEAnalysis( felems, mesh );
            if nargin < 3 || isempty(isConst)
                isConst = false;
            end
            obj.isConst = isConst;
            obj.rotations = [];
            obj.selfLoadFactor = -1;
            obj.clearWeightedStiffnessCache();
        end

        function [qfem, K] = solve(obj, x, retainStiffness)
            % solve()                  unweighted
            % solve(x)                 density-weighted
            % solve(x, retainStiffness) keeps K for a following adjoint solve
            if nargin < 2, x = []; end
            if nargin < 3 || isempty(retainStiffness), retainStiffness = false; end

            [I,J,~,~] = obj.globalMatrixIndices();
            obj.prepareRHSVectors();

            if obj.selfLoadFactor > 0
                obj.Pnodal(:) = 0;
                obj.loadElementsSelfWeight(x, 0.1);
                obj.Pfem = obj.selfLoadFactor .* obj.toFEMVector(obj.Pnodal);
                obj.Pnodal(:) = 0;
            end

            if size(obj.rotations,1) == 0
                solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
            else
                solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports), obj.rotations);
            end

            fname = obj.weightedStiffnessFunction();
            if isempty(x)
                Kvals = obj.globalMatrixAggregation('computeStifnessMatrix');
            else
                Kvals = obj.globalMatrixAggregationWeighted(fname, x);
            end

            if retainStiffness
                obj.cachedWeightedX = x(:);
                obj.cachedWeightedKvals = Kvals;
                obj.cachedWeightedFunction = fname;
            else
                obj.clearWeightedStiffnessCache();
            end

            obj.qfem = solver.solve(Kvals, obj.Pfem);
            qfem = obj.qfem;
            obj.qnodal = obj.fromFEMVector(qfem(:,1));
            if nargout > 1
                K = Kvals;      % only assembled into a matrix when asked for
            end
        end

        function lambda = solveAdjointWithLoad(obj, xPenal, P_adj_fem)
            % K(xPenal)*lambda = P_adj_fem. K is symmetric, so the adjoint
            % system is the forward system with a different RHS.
            [I, J, ~, ~] = obj.globalMatrixIndices();
            if size(obj.rotations, 1) == 0
                solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
            else
                solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports), obj.rotations);
            end
            fname = obj.weightedStiffnessFunction();
            if obj.hasCachedWeightedStiffness(xPenal, fname)
                Kvals = obj.cachedWeightedKvals;
            else
                Kvals = obj.globalMatrixAggregationWeighted(fname, xPenal);
            end
            lambda = solver.solve(Kvals, P_adj_fem);
            obj.clearWeightedStiffnessCache();
        end

        function fname = weightedStiffnessFunction(obj)
            if obj.isConst
                fname = 'computeStifnessMatrixConst';
            else
                fname = 'computeStifnessMatrix';
            end
        end

        function tf = hasCachedWeightedStiffness(obj, x, fname)
            tf = ~isempty(obj.cachedWeightedKvals) && ...
                isequal(obj.cachedWeightedFunction, fname) && ...
                numel(obj.cachedWeightedX) == numel(x) && ...
                isequal(obj.cachedWeightedX, x(:));
        end

        function clearWeightedStiffnessCache(obj)
            obj.cachedWeightedX = [];
            obj.cachedWeightedKvals = [];
            obj.cachedWeightedFunction = '';
        end

        function obj = saveMatrices(obj, filename)
            [I,J,~,~] = obj.globalMatrixIndices();
            K = sparse(I, J, obj.globalMatrixAggregationWeighted(obj.weightedStiffnessFunction(), 1));
            supports = obj.supports;
            obj.prepareRHSVectors();
            P = obj.Pfem;
            nodes = obj.mesh.nodes;
            elems = obj.mesh.elems;
            save(filename, "nodes", "elems", "K", "P", "supports");
        end
    end
end
