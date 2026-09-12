classdef SIMP_MMA_TopologyOptimizationElasticComplianceBase < SIMP_MMA_TopologyOptimization
    % Shared single-load SIMP compliance optimizer.
    %
    % stiffnessMode controls whether the objective uses the material
    % stiffness only ("linear") or material plus geometric stiffness
    % ("secondOrder").

    properties
        VolConstr;
        stiffnessMode;
    end

    methods
        function obj = SIMP_MMA_TopologyOptimizationElasticComplianceBase( ...
                Rmin, problem, penal, VolConstr, is_const, stiffnessMode)
            obj = obj@SIMP_MMA_TopologyOptimization(1, Rmin, problem, penal, is_const);
            tne = obj.FEAnalysis.getTotalElemsNumber();
            obj.gradConstrValues = zeros(1, tne);
            obj.V0 = tne * VolConstr;
            obj.VolConstr = VolConstr;
            obj.x(1:obj.totalFENumber, 1) = 0.5;
            obj.is_const = is_const;
            obj.stiffnessMode = string(stiffnessMode);
        end

        function resetAnalysis(obj)
        end

        function computeObjectiveFunctonWithGradient(obj, x)
            tne = obj.FEAnalysis.getTotalElemsNumber();
            c = 0;
            dc = zeros(tne, 1);
            xOnes = ones(tne, 1);
            xDesign = x(:);
            xPenal = xDesign .^ obj.penal;

            obj.qnodal = obj.FEAnalysis.fromFEMVector(obj.FEAnalysis.solveWeighted(xPenal));
            obj.FEAnalysis.computeElementResults(xPenal);

            ind = 1;
            for i = 1:numel(obj.FEAnalysis.felems)
                fe = obj.FEAnalysis.felems{i};
                if obj.is_const
                    fsName = 'computeStifnessMatrixConst';
                else
                    fsName = 'computeStifnessMatrix';
                end

                nelems = size(fe.elems, 1);
                nnodes = size(fe.elems, 2);
                ndofs = size(fe.eDofs, 2);
                dim = nnodes * ndofs;

                K = reshape(fe.(fsName)(obj.FEAnalysis.mesh.nodes, xOnes), dim, dim, nelems);
                if obj.stiffnessMode == "secondOrder"
                    Kg = reshape(fe.computeGeometricStifnessMatrix(obj.FEAnalysis.mesh.nodes, xOnes), ...
                        dim, dim, nelems);
                    K = K + Kg;
                end

                qelems = fe.createElemSolutionVectors(obj.qnodal);
                for j = 1:nelems
                    elemEnergy = qelems(:, j)' * K(:, :, j) * qelems(:, j);
                    c = c + xDesign(ind)^obj.penal * elemEnergy;
                    dc(ind) = -obj.penal * xDesign(ind)^(obj.penal - 1) * elemEnergy;
                    ind = ind + 1;
                end
            end

            obj.FobjValue = c;
            obj.gradFobjValue = dc;
        end

        function computeConstraintsAndGradient(obj, x)
            obj.constrValues = sum(x(:)) / obj.V0 - 1;
            obj.gradConstrValues(1, 1:obj.totalFENumber) = 1.0 / obj.V0;
        end

        function printIterationInfo(obj)
            fprintf('%5i ', obj.iteration);
            fprintf('Vrel=%5.2f ', round(obj.VolConstr * sum(obj.x) / obj.V0 * 1000) / 10);
            fprintf('constr=%6.4f ', obj.constrValues);
            fprintf('change=%6.4f ', obj.change);
            fprintf('\n');
        end
    end
end
