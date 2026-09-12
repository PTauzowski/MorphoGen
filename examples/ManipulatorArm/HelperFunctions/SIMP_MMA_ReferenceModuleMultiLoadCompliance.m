classdef SIMP_MMA_ReferenceModuleMultiLoadCompliance < SIMP_MMA_TopologyOptimization
    % Multi-load SIMP compliance optimizer for ReferenceModuleSolidModel.
    %
    % The objective is a p-norm aggregation of per-load compliances:
    %   J = (sum_k w_k * Ck^p)^(1/p)
    % or, when normalization is enabled:
    %   J = (sum_k w_k * (Ck/Ck0)^p)^(1/p)

    properties
        VolConstr;
        ReferenceModel;
        loadCases;
        loadWeights;
        pAgg;
        useComplianceNormalization;
        complianceNormalization;
        complianceValues;
        normalizedComplianceValues;
        complianceGradientValues;
        useParallel;
    end

    methods
        function obj = SIMP_MMA_ReferenceModuleMultiLoadCompliance( ...
                Rmin, referenceModel, loadCases, weights, pAgg, penal, VolConstr, is_const, useComplianceNormalization)

            obj = obj@SIMP_MMA_TopologyOptimization(1, Rmin, referenceModel.analysis, penal, is_const);
            if nargin < 9
                useComplianceNormalization = true;
            end
            if isempty(weights)
                weights = ones(numel(loadCases), 1) / numel(loadCases);
            end

            obj.ReferenceModel = referenceModel;
            obj.loadCases = loadCases(:);
            obj.loadWeights = weights(:) / sum(weights(:));
            obj.pAgg = pAgg;
            obj.VolConstr = VolConstr;
            obj.V0 = obj.totalFENumber * VolConstr;
            obj.x(1:obj.totalFENumber, 1) = VolConstr;
            obj.gradConstrValues = zeros(1, obj.totalFENumber);
            obj.useComplianceNormalization = useComplianceNormalization;
            obj.complianceNormalization = [];
            obj.complianceValues = zeros(numel(obj.loadCases), 1);
            obj.normalizedComplianceValues = zeros(numel(obj.loadCases), 1);
            obj.complianceGradientValues = zeros(obj.totalFENumber, numel(obj.loadCases));
            obj.useParallel = false;
        end

        function resetAnalysis(obj)
        end

        function initializeComplianceNormalization(obj, x)
            oldNormalization = obj.complianceNormalization;
            oldUseNormalization = obj.useComplianceNormalization;
            obj.complianceNormalization = [];
            obj.useComplianceNormalization = false;
            obj.computeObjectiveFunctonWithGradient(x);
            obj.complianceNormalization = max(obj.complianceValues, eps);
            obj.useComplianceNormalization = oldUseNormalization;
            if ~oldUseNormalization
                obj.complianceNormalization = oldNormalization;
            end
        end

        function computeObjectiveFunctonWithGradient(obj, x)
            nLoads = numel(obj.loadCases);
            C = zeros(nLoads, 1);
            dC = zeros(obj.totalFENumber, nLoads);
            qnodal = cell(nLoads, 1);
            loadPnodal = cell(nLoads, 1);

            for k = 1:nLoads
                obj.ReferenceModel.applyLoadCase(obj.loadCases{k});
                loadPnodal{k} = obj.FEAnalysis.Pnodal;
            end

            useParallel = obj.useParallel && license('test', 'Distrib_Computing_Toolbox');
            if useParallel
                FEAnalysis = obj.FEAnalysis;
                penal = obj.penal;
                is_const = obj.is_const;
                totalFENumber = obj.totalFENumber;
                parfor k = 1:nLoads
                    [C_k, dC_k, qnodal_k] = SIMP_MMA_ReferenceModuleMultiLoadCompliance.computeLoadComplianceWithGradient( ...
                        FEAnalysis, loadPnodal{k}, x, penal, is_const, totalFENumber);
                    C(k) = C_k;
                    dC(:, k) = dC_k;
                    qnodal{k} = qnodal_k;
                end
            else
                for k = 1:nLoads
                    [C(k), dC(:, k), qnodal{k}] = SIMP_MMA_ReferenceModuleMultiLoadCompliance.computeLoadComplianceWithGradient( ...
                        obj.FEAnalysis, loadPnodal{k}, x, obj.penal, obj.is_const, obj.totalFENumber);
                end
            end

            if obj.useComplianceNormalization
                if isempty(obj.complianceNormalization)
                    obj.complianceNormalization = max(C, eps);
                end
                Cagg = C ./ obj.complianceNormalization;
                dCagg = bsxfun(@rdivide, dC, obj.complianceNormalization(:)');
            else
                Cagg = C;
                dCagg = dC;
            end

            weightedSum = sum(obj.loadWeights .* (Cagg .^ obj.pAgg));
            J = weightedSum ^ (1 / obj.pAgg);
            dJ = zeros(obj.totalFENumber, 1);
            for k = 1:nLoads
                dJ = dJ + obj.loadWeights(k) * Cagg(k)^(obj.pAgg - 1) * dCagg(:, k);
            end
            dJ = J^(1 - obj.pAgg) * dJ;

            obj.complianceValues = C;
            obj.normalizedComplianceValues = Cagg;
            obj.complianceGradientValues = dC;
            obj.FobjValue = J;
            obj.gradFobjValue = dJ;
            obj.qnodal = qnodal{end};
            obj.FEAnalysis.qnodal = qnodal{end};
            obj.FEAnalysis.Pnodal = loadPnodal{end};
        end

        function computeConstraintsAndGradient(obj, x)
            obj.constrValues = sum(x(:)) / obj.V0 - 1;
            obj.gradConstrValues(1, 1:obj.totalFENumber) = 1.0 / obj.V0;
        end

        function printIterationInfo(obj)
            fprintf('%5i ', obj.iteration);
            fprintf('J=%10.4e ', obj.FobjValue);
            fprintf('Vrel=%5.2f ', round(obj.VolConstr * sum(obj.x) / obj.V0 * 1000) / 10);
            fprintf('constr=%7.4f ', obj.constrValues);
            fprintf('change=%7.4f ', obj.change);
            fprintf('\n');
        end
    end

    methods (Static, Access = private)
        function [C, dC, qnodal] = computeLoadComplianceWithGradient(FEAnalysis, Pnodal, x, penal, is_const, tne)
            C = 0;
            dC = zeros(tne, 1);
            xOnes = ones(tne, 1);
            xPhys = x(:) .^ penal;

            if is_const
                fsName = 'computeStifnessMatrixConst';
            else
                fsName = 'computeStifnessMatrix';
            end

            FEAnalysis.Pnodal = Pnodal;
            qnodal = FEAnalysis.fromFEMVector(FEAnalysis.solveWeighted(xPhys));
            FEAnalysis.computeElementResults(xPhys);

            ind = 1;
            for i = 1:numel(FEAnalysis.felems)
                fe = FEAnalysis.felems{i};
                nelems = size(fe.elems, 1);
                nnodes = size(fe.elems, 2);
                ndofs = size(fe.ndofs, 2);
                dim = nnodes * ndofs;
                K = reshape(fe.(fsName)(FEAnalysis.mesh.nodes, xOnes), dim, dim, nelems);
                qelems = fe.createElemSolutionVectors(qnodal);
                for j = 1:nelems
                    elemEnergy = qelems(:, j)' * K(:, :, j) * qelems(:, j);
                    C = C + x(ind)^penal * elemEnergy;
                    dC(ind) = -penal * x(ind)^(penal - 1) * elemEnergy;
                    ind = ind + 1;
                end
            end
        end
    end
end
