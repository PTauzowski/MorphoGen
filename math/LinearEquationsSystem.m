classdef LinearEquationsSystem < handle
        
    properties
        I,J,supports,freedofs,supdofs,newdofs;
        iK11,iK12,iK21,iK22,dim;

        Kbase_free          % sparse (nfree x nfree)
        decBase             % decomposition object (fast solves)
        baseReady (1,1) logical = false
    end
          
    methods
        
        function obj = LinearEquationsSystem(I,J,supports)
           obj.I=I;
           obj.J=J;
           obj.supports = supports;
           obj.freedofs = find(not(obj.supports));
           obj.supdofs = find(obj.supports);
           ndofs = size(obj.freedofs,1);
           nsup = size(obj.supdofs,1);
           obj.dim=ndofs+nsup;
           obj.newdofs = zeros(ndofs+nsup,1);
           obj.newdofs( obj.freedofs ) = (1:ndofs)';
           obj.newdofs( obj.supdofs ) = (1:nsup)';
           nI = obj.supports( obj.I );
           nJ = obj.supports( obj.J );
           obj.iK11 = not(nI) & not(nJ);
           obj.iK12 = not(nI) & nJ;
           obj.iK21 = nI & not(nJ);
           obj.iK22 = nI & nJ;
        end

        function q = solveClassical(obj,K,P)
             K11=sparse(obj.I,obj.J,K);
             maxK = max(max(K11));
             for k=1:max(size(obj.supdofs,1))
                K11( obj.supdofs(k), obj.supdofs(k) ) = maxK * 1.0e10;
             end
             q = K11 \ P;
        end
        function Ks = createSparseMatrix(obj,K)
            dimfree  = size( obj.freedofs, 1 );
            Ks = sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree);
        end
        function q = solve(obj,K,P)
             q=P;
             q(:)=0;
             dimfree  = size( obj.freedofs, 1 );
             q( obj.freedofs,:) = sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree) \ P( obj.freedofs,:);
        end

        function [qforms, lambdas] = solveEigenproblem(obj,K,Kg,num_eigenvalues)             
             dimfree  = size( obj.freedofs, 1 );
             [q, l] = eigs( obj.createSparseMatrix(K), obj.createSparseMatrix(Kg), num_eigenvalues, 'smallestabs');
             qforms=zeros(obj.dim,num_eigenvalues);
             qforms(obj.freedofs,:)=real(q);
             lambdas=abs(l);
        end
        
        function [q, R, error] = solveR(obj,K,P)
            q=P;
            q(:)=0;
            R=P;
            R(:)=0;
            dimfree  = size( obj.freedofs, 1 );
            dimfixed = size( obj.supdofs, 1 );
            K11=sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree);
            K21=sparse(obj.newdofs(obj.I(obj.iK21)),obj.newdofs(obj.J(obj.iK21)),K(obj.iK21),dimfixed,dimfree);
            q( obj.freedofs,:) = K11 \ P(obj.freedofs,:);
            R( obj.supdofs,:) = K21 * q( obj.freedofs,:);
            error = [ K11 * q(obj.freedofs) - P(obj.freedofs); K21 * q(obj.freedofs) - R ];
        end

        function q = solvePq(obj,K,P,q0)
             q=P;
             dimfree  = size( obj.freedofs, 1 );
             dimfixed = size( obj.supdofs, 1 );
             K11=sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree);
             K12=sparse(obj.newdofs(obj.I(obj.iK12)),obj.newdofs(obj.J(obj.iK12)),K(obj.iK12),dimfree,dimfixed);
             q( obj.freedofs,:) = K11\(P( obj.freedofs,:)-K12*q0( obj.supdofs,:));
        end

        function [q, R, error] = solvePRq(obj,K,P,q0)
            q=P;
            q(:)=0;
            R=P;
            R(:)=0;
            dimfree  = size( obj.freedofs, 1 );
            dimfixed = size( obj.supdofs, 1 );
            obj.q0 = q0( obj.supdofs,:);
            obj.P = P( obj.freedofs,:); 
            K11=sparse(obj.newdofs(obj.I(obj.iK11)),obj.newdofs(obj.J(obj.iK11)),K(obj.iK11),dimfree,dimfree);
            K12=sparse(obj.newdofs(obj.I(obj.iK12)),obj.newdofs(obj.J(obj.iK12)),K(obj.iK12),dimfree,dimfixed);
            K21=sparse(obj.newdofs(obj.I(obj.iK21)),obj.newdofs(obj.J(obj.iK21)),K(obj.iK21),dimfixed,dimfree);
            K22=sparse(obj.newdofs(obj.I(obj.iK22)),obj.newdofs(obj.J(obj.iK22)),K(obj.iK22),dimfixed,dimfixed);
            q( obj.freedofs,:) = K11\(P( obj.freedofs )-K12*q0( obj.supdofs ));
            R( obj.supdofs,:) = K21 * q( obj.freedofs ) + K22 * q0( obj.supdofs );
            error = [ K11 * q( obj.freedofs, : ) + K12*q0( obj.supdofs, : ) - P( obj.freedofs, : ); ...
                      K21 * q( obj.freedofs, : ) + K22*q0( obj.supdofs, : ) - R( obj.supdofs, : ) ];
        end

        function decBase = setBaseMatrix(obj, K)
            % Store and factorize the base free-free stiffness matrix.
            dimfree  = size(obj.freedofs, 1);
            obj.Kbase_free = sparse( ...
                obj.newdofs(obj.I(obj.iK11)), ...
                obj.newdofs(obj.J(obj.iK11)), ...
                K(obj.iK11), dimfree, dimfree);

            try
                obj.decBase = decomposition(obj.Kbase_free, 'chol');
            catch
                obj.decBase = decomposition(obj.Kbase_free, 'lu');
            end

            obj.baseReady = true;
            decBase = obj.decBase;
        end

        function freeMap = createFreeMap(obj)
            % Map global DOF index -> free index, 0 on constrained DOFs.
            freeMap = zeros(obj.dim, 1);
            freeMap(obj.freedofs) = (1:size(obj.freedofs, 1))';
        end
        
        function x = solveBase(obj, b_free, varargin)
            % Solve Kbase_free * x = b_free using cached or provided factorization.
            if nargin >= 3 && ~isempty(varargin{1})
                decBase = varargin{1};
            else
                decBase = obj.getBaseDecomposition();
            end
            x = decBase \ b_free;
        end
        
        function q = solveWoodburyFree(obj, decBase, P_free, U, C)
            % Solve q = (A + U*C*U')\P_free with Sherman-Morrison-Woodbury.
            %
            % New API:
            %   q = solveWoodburyFree(decBase, P_free, U, C)
            %
            % Backward-compatible API:
            %   q = solveWoodburyFree(P_free, U, C)
            if nargin == 4
                C = U;
                U = P_free;
                P_free = decBase;
                decBase = obj.getBaseDecomposition();
            elseif nargin == 5
                if isempty(decBase)
                    decBase = obj.getBaseDecomposition();
                end
            else
                error('solveWoodburyFree expects 3 or 4 input arguments after obj.');
            end

            y = decBase \ P_free;
            if isempty(U)
                q = y;
                return;
            end

            if size(U,1) ~= size(P_free,1)
                error('U row count must match length(P_free).');
            end
            m = size(U,2);
            if any(size(C) ~= [m, m])
                error('C must be an (m x m) matrix, where m=size(U,2).');
            end

            AU = decBase \ U;

            % Fast path for diagonal C (common for low-rank eigendecomposition).
            if norm(C - diag(diag(C)), 'fro') <= 1.0e-14 * max(1.0, norm(C, 'fro'))
                cdiag = diag(C);
                if any(abs(cdiag) <= eps)
                    error('C diagonal contains near-zero entries; cannot form C^{-1}.');
                end
                Cinv = diag(1 ./ cdiag);
            else
                Cinv = C \ eye(m);
            end

            S = Cinv + (U' * AU);
            q = y - AU * (S \ (U' * y));
        end
        
        function q = solveWoodbury(obj, decBase, P, U, C)
            % Wrapper that returns the full DOF vector with fixed DOFs zeroed.
            %
            % New API:
            %   q = solveWoodbury(decBase, P, U, C)
            %
            % Backward-compatible API:
            %   q = solveWoodbury(P, U, C)
            if nargin == 4
                C = U;
                U = P;
                P = decBase;
                decBase = obj.getBaseDecomposition();
            elseif nargin == 5
                if isempty(decBase)
                    decBase = obj.getBaseDecomposition();
                end
            else
                error('solveWoodbury expects 3 or 4 input arguments after obj.');
            end

            nFree = size(obj.freedofs, 1);
            if size(P, 1) == obj.dim
                q = zeros(size(P));
                P_free = P(obj.freedofs, :);
            elseif size(P, 1) == nFree
                q = zeros(obj.dim, size(P, 2));
                P_free = P;
            else
                error(['P must have either obj.dim rows (%d, full RHS) or ', ...
                       'nFree rows (%d, free RHS).'], obj.dim, nFree);
            end

            q(obj.freedofs, :) = obj.solveWoodburyFree(decBase, P_free, U, C);
        end

        function [q, logInfo] = solveWoodburyWithFallback( ...
                obj, decBase, K, P, U, C, info, mMax, varargin)
            % Choose Woodbury or fallback solve based on low-rank size.
            %
            % [q, logInfo] = solveWoodburyWithFallback( ...
            %   decBase, K, P, U, C, info, mMax, fallbackMode, pcgTol, pcgMaxIt)
            %
            % fallbackMode: 'direct' (default) or 'pcg'
            if nargin < 9
                error(['solveWoodburyWithFallback expects at least ', ...
                       'decBase, K, P, U, C, info and mMax.']);
            end

            if isempty(decBase)
                decBase = obj.getBaseDecomposition();
            end

            fallbackMode = 'direct';
            pcgTol = 1.0e-8;
            pcgMaxIt = 200;
            if nargin >= 10 && ~isempty(varargin{1})
                fallbackMode = lower(char(varargin{1}));
            end
            if nargin >= 11 && ~isempty(varargin{2})
                pcgTol = varargin{2};
            end
            if nargin >= 12 && ~isempty(varargin{3})
                pcgMaxIt = varargin{3};
            end

            rankU = size(U, 2);
            nChanged = NaN;
            if isstruct(info)
                if isfield(info, 'rank')
                    rankU = info.rank;
                end
                if isfield(info, 'nChanged')
                    nChanged = info.nChanged;
                end
            end

            tSolve = tic;
            pcgFlag = [];
            pcgRelRes = [];
            pcgIter = [];

            if rankU <= mMax
                q = obj.solveWoodbury(decBase, P, U, C);
                method = 'woodbury';
            else
                if isempty(K)
                    error(['K must be provided for fallback solve when ', ...
                           'info.rank exceeds mMax.']);
                end
                if strcmp(fallbackMode, 'pcg')
                    Kff = obj.createSparseMatrix(K);
                    if size(P, 1) == obj.dim
                        P_free = P(obj.freedofs, :);
                    elseif size(P, 1) == size(obj.freedofs, 1)
                        P_free = P;
                    else
                        error(['P must have either obj.dim rows (%d, full RHS) or ', ...
                               'nFree rows (%d, free RHS).'], ...
                              obj.dim, size(obj.freedofs, 1));
                    end
                    nRhs = size(P_free, 2);
                    q_free = zeros(size(P_free));
                    pcgFlag = zeros(1, nRhs);
                    pcgRelRes = zeros(1, nRhs);
                    pcgIter = zeros(1, nRhs);
                    for k = 1:nRhs
                        x0 = decBase \ P_free(:, k);
                        [qk, flagk, relresk, iterk] = pcg( ...
                            Kff, P_free(:, k), pcgTol, pcgMaxIt, [], [], x0);
                        if flagk ~= 0
                            qk = Kff \ P_free(:, k);
                        end
                        q_free(:, k) = qk;
                        pcgFlag(k) = flagk;
                        pcgRelRes(k) = relresk;
                        if numel(iterk) > 1
                            pcgIter(k) = iterk(2);
                        else
                            pcgIter(k) = iterk;
                        end
                    end
                    q = zeros(obj.dim, size(P, 2));
                    q(obj.freedofs,:) = q_free;
                    method = 'pcg';
                else
                    if size(P, 1) == obj.dim
                        Pfull = P;
                    elseif size(P, 1) == size(obj.freedofs, 1)
                        Pfull = zeros(obj.dim, size(P, 2));
                        Pfull(obj.freedofs, :) = P;
                    else
                        error(['P must have either obj.dim rows (%d, full RHS) or ', ...
                               'nFree rows (%d, free RHS).'], ...
                              obj.dim, size(obj.freedofs, 1));
                    end
                    q = obj.solve(K, Pfull);
                    method = 'direct';
                end
            end

            elapsed = toc(tSolve);
            if isnan(nChanged)
                nChangedStr = 'n/a';
            else
                nChangedStr = num2str(round(nChanged));
            end
            fprintf('[Woodbury] nChanged=%s rank=%d mMax=%d method=%s time=%.4fs\n', ...
                nChangedStr, rankU, mMax, method, elapsed);

            logInfo = struct( ...
                'rank', rankU, ...
                'nChanged', nChanged, ...
                'mMax', mMax, ...
                'method', method, ...
                'usedWoodbury', strcmp(method, 'woodbury'), ...
                'timeSec', elapsed, ...
                'pcgFlag', pcgFlag, ...
                'pcgRelRes', pcgRelRes, ...
                'pcgIter', pcgIter);
        end

        function decBase = getBaseDecomposition(obj)
            if ~obj.baseReady
                error('Base matrix not set. Call setBaseMatrix(K) first.');
            end
            decBase = obj.decBase;
        end

       
    end
end
