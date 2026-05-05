classdef StressIntensityMultiTopologyOptimization < StressIntensityTopologyOptimizationVol
    % Multi-configuration stress-intensity ESO with selectable aggregation.

    properties
        FEAnalyses;
        plLambda;
        plVol;
        plMaxStress;
        lastStableFrame;
        useParallel;
        maxdisplacement;
        intensityAggregation;
    end

    methods
        function obj = StressIntensityMultiTopologyOptimization( ...
                Rmin, FEAnalyses, maxais, penal, volFr, is_const, intensityAggregation)
            obj = obj@StressIntensityTopologyOptimizationVol(Rmin, FEAnalyses(1), ...
                maxais, penal, volFr, is_const);
            obj.FEAnalyses = FEAnalyses;
            obj.maxais = maxais;
            obj.penal = penal;
            obj.x(:) = 1;
            obj.xfat = obj.x;
            obj.pnormfat = 100;
            obj.plLambda = [];
            obj.plVol = [];
            obj.plMaxStress = [];
            obj.useParallel = false;
            obj.maxdisplacement = [];
            obj.intensityAggregation = string(intensityAggregation);
        end

        function printIterationInfo(obj)
            fprintf('%5i ', obj.iteration);
            fprintf('Vrel=%2.1f ', round(sum(obj.x) / obj.V0 * 1000) / 10);
            fprintf('lambda=%5.3g ', obj.FEAnalysis.lambda);
            fprintf('\n');
            obj.plLambda = [obj.plLambda abs(obj.FEAnalysis.lambda)];
            obj.plVol = [obj.plVol round(sum(obj.x) / obj.V0 * 1000) / 10];
            if abs(obj.FEAnalysis.lambda) >= 1
                obj.lastStableFrame = obj.iteration;
            end
        end

        function ais = computeAverageIntensities(obj)
            nFEAnalyses = numel(obj.FEAnalyses);
            nElems = obj.FEAnalysis.getTotalElemsNumber();
            sais = zeros(nElems, nFEAnalyses);
            maxstress = zeros(1, nFEAnalyses);
            maxdisplacement = zeros(1, nFEAnalyses);
            qnodal = cell(1, nFEAnalyses);

            FEAnalyses = obj.FEAnalyses;
            elem_inds = obj.elem_inds;
            xPenal = obj.x .^ obj.penal;

            useParallel = obj.useParallel && license('test', 'Distrib_Computing_Toolbox');
            if useParallel
                parfor k = 1:nFEAnalyses
                    [sais_k, maxstress_k, maxdisplacement_k, qnodal_k] = ...
                        StressIntensityMultiTopologyOptimization.computeConfigStressIntensity( ...
                        FEAnalyses(k), elem_inds, xPenal, nElems);
                    sais(:, k) = sais_k;
                    maxstress(k) = maxstress_k;
                    maxdisplacement(k) = maxdisplacement_k;
                    qnodal{k} = qnodal_k;
                end
            else
                for k = 1:nFEAnalyses
                    [sais_k, maxstress_k, maxdisplacement_k, qnodal_k] = ...
                        StressIntensityMultiTopologyOptimization.computeConfigStressIntensity( ...
                        FEAnalyses(k), elem_inds, xPenal, nElems);
                    sais(:, k) = sais_k;
                    maxstress(k) = maxstress_k;
                    maxdisplacement(k) = maxdisplacement_k;
                    qnodal{k} = qnodal_k;
                end
            end

            obj.qnodal = qnodal{end};
            obj.maxstress = [obj.maxstress maxstress];
            obj.maxdisplacement = [obj.maxdisplacement maxdisplacement];
            obj.plMaxStress = [obj.plMaxStress; maxstress];

            switch obj.intensityAggregation
                case "max"
                    ais = max(sais, [], 2);
                case {"sum", "average", "avg"}
                    ais = sum(sais, 2);
                otherwise
                    error('StressIntensityMultiTopologyOptimization:UnknownAggregation', ...
                        'Unknown stress-intensity aggregation "%s".', obj.intensityAggregation);
            end
            ais = ais / max(ais);
        end
    end

    methods (Static)
        function [sais_k, maxstress_k, maxdisplacement_k, qnodal_k] = ...
                computeConfigStressIntensity(analysis, elem_inds, xPenal, nElems)
            qfem = analysis.solveWeighted(xPenal);
            qnodal_k = analysis.fromFEMVector(qfem(:, 1));
            analysis.computeElementResults(xPenal);

            sais_k = zeros(nElems, 1);
            for i = 1:size(analysis.felems, 2)
                hmIndex = find(analysis.felems{i}.results.names == "sHM");
                for j = 1:size(analysis.felems{i}.elems, 1)
                    sais_k(elem_inds{i}(j)) = mean( ...
                        analysis.felems{i}.results.nodal.all( ...
                        analysis.felems{i}.elems(j, :), hmIndex));
                end
            end

            maxstress_k = max(sais_k);
            displacementDofs = analysis.findDOFsIndices(["ux", "uy", "uz"]);
            maxdisplacement_k = max(vecnorm(qnodal_k(:, displacementDofs), 2, 2));
        end
    end
end
