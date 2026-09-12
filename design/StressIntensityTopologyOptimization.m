classdef (Abstract) StressIntensityTopologyOptimization < TopologyOptimization
    
    properties
        maxstress, ais, elem_list, maxais, penal, xr, xfat, pnormfat, allx;
        useWoodburyReanalysis logical = false;
        woodburyEigTol = 1.0e-10;
        woodburyRankMax = 300;
        woodburyFallbackMode = 'direct'; % 'direct' or 'pcg'
        woodburyVerbose logical = false;
        woodburyMonitorPrint logical = false;
        woodburyRebaseOnFallback logical = true;
        woodburyState;
        woodburyMonitorHistory;
    end

    methods
            
        function obj = StressIntensityTopologyOptimization(numberOfConstraints,Rmin,FEAnalysis,maxais,penal,is_const)
            obj=obj@TopologyOptimization(numberOfConstraints,Rmin,FEAnalysis,is_const);
            obj.maxais=maxais;
            obj.penal=penal;
            obj.x(:)=1;
            obj.xfat=obj.x;
            obj.pnormfat=100;
            obj.resetWoodburyState();
            obj.resetWoodburyMonitor();
        end

        function ret = isNotFinished(obj)
            obj.constrValues = obj.computeInequalityConstraints(obj.x);
            if find( obj.constrValues < 0 , 1 )
                    ret = false;
                    return;
            end
            ret=true;
        end

        function resetAnalysis(obj)
            obj.x(:)=1;
            obj.totalFENumber = obj.FEAnalysis.getTotalElemsNumber();
            obj.erased_elems = false(obj.totalFENumber,1);
            obj.resetWoodburyState();
            obj.resetWoodburyMonitor();
        end

        function updateDesign(obj)
            obj.ais = obj.weights * obj.computeAverageIntensities();
            obj.removeStressed();
            wais=obj.ais/max(obj.ais);
            obj.x(obj.erased_elems) = min(1, max( 0.001, obj.x(obj.erased_elems) .* wais(obj.erased_elems) .^ obj.penal)); 
            obj.x( not(obj.erased_elems) ) = 1; 
            obj.x( obj.const_elems ) = 1;
            obj.allx=[obj.allx obj.x];
        end

        function removeStressed( obj )
            a=20;
            b=0;
            tc1=0.05;
            tc2=0.005;
            y5 = (1 - sum( obj.x )/obj.V0)/(1-obj.Vend) * 0.05 + (sum( obj.x )/obj.V0)/(1-obj.Vend) * 0.005;
            y6 = (1 - sum( obj.x )/obj.V0)/(1-obj.Vend) * 0.005 + (sum( obj.x )/obj.V0)/(1-obj.Vend) * 0.05;

            x = 1-sum( obj.x )/obj.V0;

            y1=abs(tc2-tc1)*(2./(1+exp(-(a*(x-b))))-1)+min(tc1,tc2);
            y2=abs(tc2-tc1)*(1-(2./(1+exp(-(a*(x-b))))-1))+min(tc1,tc2);

            a=20;
            b=0.2;
            y3=abs(tc2-tc1)*(1./(1+exp(-(a*(x-b)))))+min(tc1,tc2);
            y4=abs(tc2-tc1)*(1-(1./(1+exp(-(a*(x-b))))))+min(tc1,tc2);

            %obj.maxais = y4;

            max_elem_removal_factor = 0.025;

            notErasedID = find( not( obj.erased_elems )  );
            notErasedID = setxor(notErasedID,intersect(notErasedID, obj.const_elems));
            maxaisprc =(max(obj.ais(notErasedID)) - min(obj.ais(notErasedID)) ) * obj.maxais;
            obj.elem_list = obj.ais(notErasedID) < min(obj.ais(notErasedID)) + maxaisprc;
            if size(find(obj.elem_list),1)>size(obj.ais,1)*max_elem_removal_factor
                    [~, ai]=sort(obj.ais(notErasedID));
                    obj.elem_list=false(size(obj.ais(notErasedID),1),1);
                    obj.elem_list(ai(1:round(size(obj.ais,1)*max_elem_removal_factor)))=true;
            end

            obj.elem_list = notErasedID( obj.elem_list );
            obj.erased_elems( obj.elem_list ) = true;
        end

        function ais = computeAverageIntensities(obj)
            qfem = obj.solveFEMForCurrentDensity((obj.x).^obj.penal);
            n_rhsv = size(qfem,2);
            hm_stress = zeros(obj.FEAnalysis.getTotalElemsNumber(),n_rhsv);
            for k=1:n_rhsv
                obj.FEAnalysis.qnodal = obj.FEAnalysis.fromFEMVector(qfem(:,k));
                obj.FEAnalysis.computeElementResults(obj.x.^obj.penal);
                for i=1:size(obj.FEAnalysis.felems,2)
                   hmIndex=find(obj.FEAnalysis.felems{i}.results.names == "sHM");
                   for j=1:size(obj.FEAnalysis.felems{i}.elems,1)
                        %ais(obj.elem_inds{i}(j)) = mean( obj.linearElasticProblem.felems{i}.results.GPvalues(hmIndex,j,:) );
                        hm_stress(obj.elem_inds{i}(j),k) = mean( obj.FEAnalysis.felems{i}.results.nodal.all(obj.FEAnalysis.felems{i}.elems(j,:),hmIndex) );
                   end
                end
            end
            max_hmstress = max(hm_stress,[],"all");
            obj.maxstress = [ obj.maxstress max_hmstress ];
            ais = max(hm_stress,[],2) / max_hmstress;
        end

        % function ais = computeAverageIntensities(obj)
        %     obj.qnodal = obj.FEAnalysis.solve((obj.x).^obj.penal);
        %     obj.FEAnalysis.computeElementResults(obj.x.^obj.penal);
        %     ais = zeros(obj.FEAnalysis.getTotalElemsNumber(),1);
        %     for i=1:size(obj.FEAnalysis.felems,2)
        %        hmIndex=find(obj.FEAnalysis.felems{i}.results.names == "sHM");
        %        for j=1:size(obj.FEAnalysis.felems{i}.elems,1)
        %             %ais(obj.elem_inds{i}(j)) = mean( obj.linearElasticProblem.felems{i}.results.GPvalues(hmIndex,j,:) );
        %             ais(obj.elem_inds{i}(j)) = mean( obj.FEAnalysis.felems{i}.results.nodal.all(obj.FEAnalysis.felems{i}.elems(j,:),hmIndex) );
        %        end
        %     end
        %     obj.maxstress = [ obj.maxstress max(ais) ];
        %     ais = ais / max(ais);
        % end

        function setFrame( obj, k )
            if k > 0 && k <= size( obj.allx,2 )
                obj.x=obj.allx(:,k);
                obj.FEAnalysis.computeElementResults(obj.x.^obj.penal);
                obj.iteration = k;
                return;
            end
            fprintf("The specified iteration number %d is outside the allowed range %d - %d\n",k,1,size( obj.allx,2 ));
        end

        function ifound = findFrame(obj, vol )
            i1=1; 
            i2=size(obj.allx,2);
            v1 = sum(obj.allx(:,1))/size(obj.allx,1);
            v2 = sum(obj.allx(:,end))/size(obj.allx,1);
            if vol > v1 || vol < v2
                fprintf("The specified volume fraction %1.3f is outside the existing range %1.3f - %1.3f\n",vol,v1,v2);
            else
                while true
                    i=round((i1+i2)/2);
                    v=sum(obj.allx(:,i))/size(obj.allx,1);
                    %fprintf("Mid point i=%d, v=%1.3f\n",i,v);
                    if i2-i1==1
                        ifound=i1;
                        break;
                    end
                    if v>vol
                        v1=v;
                        i1=i;
                    else
                        v2=v;
                        i2=i;
                    end
                end
            end
            
        end

        function printVolumeFractionRanges(obj)
            i1=1; 
            i2=size(obj.allx,2);
            obj.setFrame(i1);
            [v1, va1, vc1] = obj.computeVolumeFraction();
            obj.setFrame(i2);
            [v2, va2, vc2] = obj.computeVolumeFraction();
            fprintf("Frame %5d, VolFr=%1.3f, VolFr_active=%1.3f, VolFr_const=%1.3f\n",i1,v1,va1,vc1);
            fprintf("Frame %5d, VolFr=%1.3f, VolFr_active=%1.3f, VolFr_const=%1.3f\n",i2,v2,va2,vc2);
        end

        function configureWoodburyReanalysis(obj, useFlag, varargin)
            % Enable/disable Woodbury reanalysis path for eligible stress topology runs.
            % Optional args:
            %   rankMax, eigTol, fallbackMode, verbose, monitorPrint
            obj.useWoodburyReanalysis = logical(useFlag);

            if nargin >= 3 && ~isempty(varargin{1})
                obj.woodburyRankMax = varargin{1};
            end
            if nargin >= 4 && ~isempty(varargin{2})
                obj.woodburyEigTol = varargin{2};
            end
            if nargin >= 5 && ~isempty(varargin{3})
                obj.woodburyFallbackMode = lower(char(varargin{3}));
            end
            if nargin >= 6 && ~isempty(varargin{4})
                obj.woodburyVerbose = logical(varargin{4});
            end
            if nargin < 7
                obj.woodburyMonitorPrint = obj.woodburyVerbose;
            elseif ~isempty(varargin{5})
                obj.woodburyMonitorPrint = logical(varargin{5});
            end

            obj.resetWoodburyState();
            obj.resetWoodburyMonitor();
        end

        function T = getWoodburyMonitorTable(obj)
            % Return per-iteration Woodbury diagnostics as a table.
            if isempty(obj.woodburyMonitorHistory)
                T = table();
                return;
            end
            T = struct2table(obj.woodburyMonitorHistory);
        end

        function printWoodburyMonitorSummary(obj, lastN)
            % Print a compact monitoring summary.
            if nargin < 2 || isempty(lastN)
                lastN = 10;
            end
            if isempty(obj.woodburyMonitorHistory)
                fprintf('[StressTO-WB] monitor: no records.\n');
                return;
            end
            T = obj.getWoodburyMonitorTable();
            n = height(T);
            i0 = max(1, n - lastN + 1);
            fprintf('[StressTO-WB] monitor summary (last %d of %d):\n', n - i0 + 1, n);
            disp(T(i0:n, :));
        end

    end

    methods (Access = protected)
        function qfem = solveFEMForCurrentDensity(obj, wCurr)
            tSolve = tic;
            usedWoodbury = false;
            method = 'direct';
            rankVal = 0;
            nChanged = 0;

            if ~obj.useWoodburyReanalysis
                qfem = obj.FEAnalysis.solve(wCurr);
                return;
            end

            if obj.woodburyState.eligibilityChecked && ~obj.woodburyState.eligible
                qfem = obj.FEAnalysis.solve(wCurr);
                method = 'direct-not-eligible';
                if obj.woodburyVerbose && ~obj.woodburyMonitorPrint && ~obj.woodburyState.notEligiblePrinted
                    fprintf(['[StressTO-WB] disabled for this run: eligibility check failed ', ...
                             '(requires single element block with fixed local DOFs, e.g. Q4/H8).\n']);
                    obj.woodburyState.notEligiblePrinted = true;
                end
                obj.recordWoodburyMonitor(method, usedWoodbury, rankVal, nChanged, toc(tSolve));
                return;
            end

            try
                [qfem, usedWoodbury, info] = obj.tryWoodburyReanalysisSolve(wCurr);
                if isfield(info, 'method')
                    method = info.method;
                end
                if isfield(info, 'rank')
                    rankVal = info.rank;
                end
                if isfield(info, 'nChanged')
                    nChanged = info.nChanged;
                end

                if obj.woodburyVerbose && ~obj.woodburyMonitorPrint
                    if usedWoodbury
                        fprintf('[StressTO-WB] used Woodbury rank=%d nChanged=%d\n', ...
                            rankVal, nChanged);
                    else
                        if strcmp(method, 'direct-not-eligible')
                            if ~obj.woodburyState.notEligiblePrinted
                                fprintf(['[StressTO-WB] disabled for this run: eligibility check failed ', ...
                                         '(requires single element block with fixed local DOFs, e.g. Q4/H8).\n']);
                                obj.woodburyState.notEligiblePrinted = true;
                            end
                        else
                            fprintf('[StressTO-WB] fallback method=%s rank=%d nChanged=%d\n', ...
                                method, rankVal, nChanged);
                        end
                    end
                end
            catch ME
                if isempty(ME.stack)
                    warning('StressIntensityTopologyOptimization:WoodburyFallback', ...
                        'Woodbury path failed (%s). Falling back to direct solve.', ME.message);
                else
                    warning('StressIntensityTopologyOptimization:WoodburyFallback', ...
                        'Woodbury path failed (%s). Falling back to direct solve. [%s:%d]', ...
                        ME.message, ME.stack(1).name, ME.stack(1).line);
                end
                obj.resetWoodburyState();
                qfem = obj.FEAnalysis.solve(wCurr);
                method = 'direct-error-fallback';
            end

            obj.recordWoodburyMonitor(method, usedWoodbury, rankVal, nChanged, toc(tSolve));
        end

        function [qfem, usedWoodbury, infoOut] = tryWoodburyReanalysisSolve(obj, wCurr)
            wCurr = wCurr(:);

            if ~obj.isWoodburyEligible(wCurr)
                obj.woodburyState.eligibilityChecked = true;
                obj.woodburyState.eligible = false;
                qfem = obj.FEAnalysis.solve(wCurr);
                usedWoodbury = false;
                infoOut = struct('rank', 0, 'nChanged', 0, 'method', 'direct-not-eligible');
                return;
            end
            obj.woodburyState.eligibilityChecked = true;
            obj.woodburyState.eligible = true;

            if ~obj.woodburyState.initialized
                [I, J, ~] = obj.FEAnalysis.globalMatrixIndices();
                solver = LinearEquationsSystem(I, J, obj.FEAnalysis.toFEMVector(obj.FEAnalysis.supports));
                Kbase = obj.FEAnalysis.assemblyGlobalMatrix('computeStifnessMatrix', wCurr, obj.is_const);
                decBase = solver.setBaseMatrix(Kbase);

                obj.woodburyState.initialized = true;
                obj.woodburyState.solver = solver;
                obj.woodburyState.decBase = decBase;
                obj.woodburyState.baseW = wCurr;
                obj.woodburyState.freeMap = solver.createFreeMap();
                obj.woodburyState.nFree = numel(solver.freedofs);
            end

            state = obj.woodburyState;
            Pcurr = obj.getCurrentRightHandSideForWoodbury(wCurr);
            changeTol = eps(max(1.0, norm(wCurr, inf)));
            changedGlobal = find(abs(wCurr - state.baseW) > changeTol);

            if isempty(changedGlobal)
                qfem = state.solver.solveWoodbury(state.decBase, Pcurr, ...
                    zeros(state.nFree, 0), zeros(0, 0));
                usedWoodbury = true;
                infoOut = struct('rank', 0, 'nChanged', 0, 'method', 'woodbury');
                return;
            end

            [changedUsed, changedLocal] = obj.mapChangedToLocalElements(changedGlobal);
            if isempty(changedLocal)
                qfem = obj.FEAnalysis.solve(wCurr);
                usedWoodbury = false;
                infoOut = struct('rank', 0, 'nChanged', 0, 'method', 'direct-no-local-map');
                return;
            end

            edofs_list = obj.buildElementDofsList(changedLocal);
            Ke0_vec_list = obj.buildElementKeVecList(changedLocal);
            [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
                edofs_list, Ke0_vec_list, state.baseW(changedUsed), wCurr(changedUsed), ...
                state.freeMap, state.nFree, obj.woodburyEigTol);

            if info.rank <= obj.woodburyRankMax
                qfem = state.solver.solveWoodbury(state.decBase, Pcurr, U, C);
                usedWoodbury = true;
                infoOut = info;
                infoOut.method = 'woodbury';
            else
                Kcurr = obj.FEAnalysis.assemblyGlobalMatrix('computeStifnessMatrix', wCurr, obj.is_const);
                if strcmpi(obj.woodburyFallbackMode, 'pcg')
                    [qfem, logInfo] = state.solver.solveWoodburyWithFallback( ...
                        state.decBase, Kcurr, Pcurr, U, C, info, ...
                        obj.woodburyRankMax, 'pcg');
                    infoOut = info;
                    infoOut.method = logInfo.method;
                else
                    qfem = state.solver.solve(Kcurr, Pcurr);
                    infoOut = info;
                    infoOut.method = 'direct';
                end

                usedWoodbury = false;
                if obj.woodburyRebaseOnFallback
                    state.decBase = state.solver.setBaseMatrix(Kcurr);
                    state.baseW = wCurr;
                    obj.woodburyState = state;
                end
            end
        end

        function ok = isWoodburyEligible(obj, wCurr)
            ok = false;
            if ~exist('lowRankFromElementUpdates_vecKe_array', 'file')
                return;
            end
            if size(obj.FEAnalysis.felems, 2) ~= 1
                return;
            end
            if numel(obj.elem_inds) ~= 1
                return;
            end
            if isempty(wCurr)
                return;
            end
            fe = obj.FEAnalysis.felems{1};
            if numel(obj.elem_inds{1}) ~= numel(wCurr)
                return;
            end
            nNodesPerElem = size(fe.elems, 2);
            nDofPerNode = size(obj.FEAnalysis.ndofs, 2);
            if nNodesPerElem * nDofPerNode <= 0
                return;
            end
            if numel(wCurr) ~= size(fe.elems, 1)
                return;
            end
            ok = true;
        end

        function edofs_list = buildElementDofsList(obj, changedLocal)
            fe = obj.FEAnalysis.felems{1};
            nNodesPerElem = size(fe.elems, 2);
            nElemDofs = size(fe.eDofs, 2);
            nGlobalDofs = size(obj.FEAnalysis.ndofs, 2);
            [~, ~, idofs] = intersect(fe.eDofs, obj.FEAnalysis.ndofs);
            if numel(idofs) ~= nElemDofs
                error('Woodbury reanalysis: could not map element DOFs to global DOF ordering.');
            end
            allEdofs = reshape( ...
                (repelem(fe.elems, 1, nElemDofs) - 1) * nGlobalDofs + ...
                repmat(idofs', size(fe.elems, 1), nNodesPerElem), ...
                size(fe.elems, 1), nNodesPerElem * nElemDofs);
            edofs_list = allEdofs(changedLocal, :);
            if isempty(edofs_list) || size(edofs_list, 2) == 0
                error('Woodbury reanalysis: invalid local DOF connectivity.');
            end
        end

        function Ke0_vec_list = buildElementKeVecList(obj, changedLocal)
            fe = obj.FEAnalysis.felems{1};
            Ke0 = fe.computeStifnessMatrix(obj.FEAnalysis.mesh.nodes, changedLocal(:));
            if ndims(Ke0) == 2
                Ke0 = reshape(Ke0, size(Ke0,1), size(Ke0,2), 1);
            end
            Ke0_vec_list = reshape(permute(Ke0, [3, 1, 2]), numel(changedLocal), []);
            if size(Ke0, 1) ~= size(Ke0, 2)
                error('Woodbury reanalysis expects square local element stiffness matrices.');
            end
            if size(Ke0_vec_list, 2) ~= size(Ke0, 1)^2
                error('Woodbury reanalysis expects baseline element matrices reshaped to ndLocal^2 entries.');
            end
        end

        function [changedUsed, changedLocal] = mapChangedToLocalElements(obj, changedGlobal)
            elemGlobal = obj.elem_inds{1}(:);
            [tf, loc] = ismember(changedGlobal(:), elemGlobal);
            changedUsed = changedGlobal(tf);
            changedLocal = loc(tf);
        end

        function P = getCurrentRightHandSideForWoodbury(obj, wCurr)
            % Keep RHS update consistent with FEAnalysis.solve(...) for self-weight runs.
            if isprop(obj.FEAnalysis, 'selfLoadFactor') && ~isempty(obj.FEAnalysis.selfLoadFactor) ...
                    && obj.FEAnalysis.selfLoadFactor > 0
                obj.FEAnalysis.Pnodal(:) = 0;
                obj.FEAnalysis.loadElementsSelfWeight(wCurr, 0.1);
                P = obj.FEAnalysis.selfLoadFactor .* obj.FEAnalysis.toFEMVector(obj.FEAnalysis.Pnodal);
                obj.FEAnalysis.Pfem = P;
                obj.FEAnalysis.Pnodal(:) = 0;
            else
                P = obj.FEAnalysis.Pfem;
            end

            if isempty(P)
                error(['Woodbury reanalysis: RHS vector is empty. ', ...
                    'Set FEAnalysis.Pfem or enable selfLoadFactor.']);
            end
        end

        function resetWoodburyState(obj)
            obj.woodburyState = struct( ...
                'initialized', false, ...
                'solver', [], ...
                'decBase', [], ...
                'baseW', [], ...
                'freeMap', [], ...
                'nFree', 0, ...
                'eligibilityChecked', false, ...
                'eligible', false, ...
                'notEligiblePrinted', false);
        end

        function resetWoodburyMonitor(obj)
            obj.woodburyMonitorHistory = struct( ...
                'iteration', {}, ...
                'method', {}, ...
                'usedWoodbury', {}, ...
                'rank', {}, ...
                'nChanged', {}, ...
                'timeSec', {}, ...
                'vrel', {}, ...
                'eligible', {});
        end

        function recordWoodburyMonitor(obj, method, usedWoodbury, rankVal, nChanged, timeSec)
            rec = struct( ...
                'iteration', obj.iteration, ...
                'method', char(method), ...
                'usedWoodbury', logical(usedWoodbury), ...
                'rank', rankVal, ...
                'nChanged', nChanged, ...
                'timeSec', timeSec, ...
                'vrel', sum(obj.x) / max(1, obj.V0), ...
                'eligible', obj.woodburyState.eligible);
            obj.woodburyMonitorHistory(end+1,1) = rec;

            if obj.woodburyMonitorPrint
                if strcmp(rec.method, 'direct-not-eligible')
                    if ~obj.woodburyState.notEligiblePrinted
                        fprintf(['[StressTO-WB] disabled for this run: eligibility check failed ', ...
                                 '(requires single element block with fixed local DOFs, e.g. Q4/H8).\n']);
                        obj.woodburyState.notEligiblePrinted = true;
                    end
                else
                    fprintf('[StressTO-WB] it=%d method=%s rank=%d nChanged=%d t=%.4fs\n', ...
                        rec.iteration, rec.method, rec.rank, rec.nChanged, rec.timeSec);
                end
            end
        end
    end
end
