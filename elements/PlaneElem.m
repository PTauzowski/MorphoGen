classdef PlaneElem < FiniteElement
      
    methods
        function obj = PlaneElem(sf, elems)
            obj = obj@FiniteElement(sf,elems);
        end
        function obj = setThickness( obj, h )
            obj.props.h = h;
        end
        function N = shapeMatrix( obj, points )
            nnodes = size(obj.elems, 2 );
            nd = size( obj.ndofs, 2 );
            np = size(points,1);
            N = zeros( nd, nd*nnodes, np );
            Nsf = obj.sf.computeValue(points);
            for k=1:np
                N(1,1:2:nd*nnodes-1,k) = Nsf(k,:);
                N(2,2:2:nd*nnodes,k) = Nsf(k,:);
            end
        end
        function [Jinv, detJ] = jacobianInversePages(obj, nodes, elemIds, dNtr)
            elemX = obj.elementNodePages(nodes, elemIds);
            nelems = numel(elemIds);
            nip = size(dNtr, 3);
            detJ = zeros(1, 1, nelems, nip);
            Jinv = zeros(2, 2, nelems, nip);
            for ip = 1:nip
                J = pagemtimes(dNtr(:, :, ip), elemX);
                detJ(:, :, :, ip) = J(1,1,:) .* J(2,2,:) - J(1,2,:) .* J(2,1,:);
                dj = detJ(:, :, :, ip);
                Jinv(1,1,:,ip) =  J(2,2,:) ./ dj;
                Jinv(1,2,:,ip) = -J(1,2,:) ./ dj;
                Jinv(2,1,:,ip) = -J(2,1,:) ./ dj;
                Jinv(2,2,:,ip) =  J(1,1,:) ./ dj;
            end
        end
        function B = strainBPages(obj, Jinv, dNtr)
            nnodes = size(obj.elems, 2);
            dim = 2 * nnodes;
            nelems = size(Jinv, 3);
            nip = size(Jinv, 4);
            B = zeros(3, dim, nelems, nip);
            for ip = 1:nip
                dNx = pagemtimes(Jinv(:, :, :, ip), dNtr(:, :, ip));
                cols = 1:nnodes;
                c1 = 2 * cols - 1;
                c2 = 2 * cols;
                B(1, c1, :, ip) = dNx(1, cols, :);
                B(2, c2, :, ip) = dNx(2, cols, :);
                B(3, c1, :, ip) = dNx(2, cols, :);
                B(3, c2, :, ip) = dNx(1, cols, :);
            end
        end
        function Sg = geometricGradientPages(~, Jinv, dNtr)
            nelems = size(Jinv, 3);
            nip = size(Jinv, 4);
            nnodes = size(dNtr, 2);
            Sg = zeros(2, nnodes, nelems, nip);
            for ip = 1:nip
                Sg(:, :, :, ip) = pagemtimes(Jinv(:, :, :, ip), dNtr(:, :, ip));
            end
        end
        function s = geometricStressPages(obj, elemIds, nip)
            nelems = numel(elemIds);
            s = zeros(2, 2, nelems, nip);
            stress = obj.results.gp.stress(:, elemIds, :);
            s(1,1,:,:) = reshape(stress(1,:,:), 1, 1, nelems, nip);
            s(2,2,:,:) = reshape(stress(2,:,:), 1, 1, nelems, nip);
            s(1,2,:,:) = reshape(stress(3,:,:), 1, 1, nelems, nip);
            s(2,1,:,:) = reshape(stress(3,:,:), 1, 1, nelems, nip);
        end
        function K = composeGeometricPages(~, So, scale)
            nnodes = size(So, 1);
            nelems = size(So, 3);
            dim = 2 * nnodes;
            K = zeros(dim, dim, nelems);
            So = So .* reshape(scale, 1, 1, []);
            K(1:2:dim, 1:2:dim, :) = So;
            K(2:2:dim, 2:2:dim, :) = So;
        end
        function K = computeStifnessMatrix(obj, nodes, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            x = obj.elementScale(nelems, varargin{:});
            K = zeros(dim, dim, nelems);
            chunkSize = obj.assemblyChunkSize(dim);
            D = obj.mat.D;
            h = obj.props.h;
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [Jinv, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                B = obj.strainBPages(Jinv, dNtr);
                K(:, :, elemIds) = obj.integratePagematrix(B, D, detJ, ...
                    integrator.weights, h * x(elemIds));
            end
            K = obj.flattenElementMatrices(K);
        end
        function K = computeGeometricStifnessMatrix(obj, nodes, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            x = obj.elementScale(nelems, varargin{:});
            K = zeros(dim, dim, nelems);
            chunkSize = obj.assemblyChunkSize(dim);
            h = obj.props.h;
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [Jinv, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                Sg = obj.geometricGradientPages(Jinv, dNtr);
                s = obj.geometricStressPages(elemIds, size(integrator.points, 1));
                So = obj.integratePagematrix(Sg, s, detJ, integrator.weights, ones(numel(elemIds), 1));
                K(:, :, elemIds) = obj.composeGeometricPages(So, h * x(elemIds));
            end
            K = obj.flattenElementMatrices(K);
        end
        function M = computeMassMatrix(obj, nodes, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            N = obj.shapeMatrix(integrator.points);
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            x = obj.elementScale(nelems, varargin{:});
            M = zeros(dim, dim, nelems);
            chunkSize = obj.assemblyChunkSize(dim);
            h = obj.props.h;
            if isprop(obj.mat, 'M') && ~isempty(obj.mat.M)
                massMatrix = obj.mat.M;
            elseif isprop(obj.mat, 'rho') && ~isempty(obj.mat.rho)
                massMatrix = obj.mat.rho * eye(ndofs);
            else
                massMatrix = eye(ndofs);
            end
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [~, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                nc = numel(elemIds);
                Npages = zeros(ndofs, dim, nc, size(integrator.points, 1));
                for ip = 1:size(integrator.points, 1)
                    Npages(:, :, :, ip) = repmat(N(:, :, ip), 1, 1, nc);
                end
                M(:, :, elemIds) = obj.integratePagematrix(Npages, massMatrix, detJ, ...
                    integrator.weights, h * x(elemIds));
            end
            M = obj.flattenElementMatrices(M);
        end
        function dK = computeStifnessMatrixGradMat(obj, nodes, q, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            nsens = size(obj.mat.dD,3);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            x = obj.elementScale(nelems, varargin{:});
            dK = zeros( dim, nelems, nsens );
            h = obj.props.h;
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            chunkSize = obj.assemblyChunkSize(dim);
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [Jinv, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                B = obj.strainBPages(Jinv, dNtr);
                qpages = reshape(qelems(:, elemIds), dim, 1, []);
                for s=1:nsens
                    Ke = obj.integratePagematrix(B, obj.mat.dD(:,:,s), detJ, ...
                        integrator.weights, h * x(elemIds));
                    dK(:, elemIds, s) = squeeze(pagemtimes(Ke, qpages));
                end
            end
            dK=reshape(dK,[nnodes * ndofs*nelems nsens]);
        end
        function dK = computeStifnessMatrixGradX(obj, nodes, q, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            x = obj.elementScale(nelems, varargin{:});
            dK = zeros( dim, nelems );
            h = obj.props.h;
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            D=obj.mat.D;
            chunkSize = obj.assemblyChunkSize(dim);
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [Jinv, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                B = obj.strainBPages(Jinv, dNtr);
                Ke = obj.integratePagematrix(B, D, detJ, integrator.weights, h * x(elemIds));
                qpages = reshape(qelems(:, elemIds), dim, 1, []);
                dK(:, elemIds) = squeeze(pagemtimes(Ke, qpages));
            end
        end
        function K = computeStifnessMatrixConst(obj, nodes, x)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            D = obj.mat.D;
            h = obj.props.h;
            x = obj.elementScale(nelems, x);
            [Jinv, detJ] = obj.jacobianInversePages(nodes, 1, dNtr);
            B = obj.strainBPages(Jinv, dNtr);
            Ke = obj.integratePagematrix(B, D, detJ, integrator.weights, h);
            K = Ke .* reshape(x, 1, 1, []);
            K = obj.flattenElementMatrices(K);
        end
        function [P, volume] = loadLineIntegral(obj, mode, nodes, edges, dofnames, di, P, valueFn)
            inds = obj.findDofIndices( dofnames );
            integrator = obj.sf.edgesf.createIntegrator();
            N = obj.sf.edgesf.computeValue( integrator.points );
            dN = obj.sf.edgesf.computeGradient( integrator.points );
            nip = size(dN,2);
            np  = size(dN,1);
            for k=1:size(edges,1)
                elemX = nodes(edges(k,:),:);
                dXY = dN * elemX;
                Pfn = zeros( size(elemX,1), size( obj.ndofs,2) );
                Pfn(:,inds) = valueFn( N*elemX );
                if mode =="local"
                    nXY = dXY./vecnorm(dXY')';
                    tXY = [ nXY(:,2) -nXY(:,1) ];
                    Pfn = nXY .* Pfn(:,1) + tXY .* Pfn(:,2);
                end
                Pe = zeros( size(Pfn) );
                volume = 0;
                for i=1:nip
                    dt = sum( dXY(i,:) .* dXY(i,:) );
                    Pe = Pe + integrator.weights( i ) * sqrt( dt ) * N(i,:)' .* Pfn(i,:);
                    volume = volume + integrator.weights( i ) * sqrt( dt );
                end
                P(edges(k,:),di) = P(edges(k,:),di) + Pe(:,inds);
            end

        end
        
        function plotWired(obj, nodes, varargin)
            hold on;
            daspect([1 1 1]);
            if nargin==2
                patch('Vertices', nodes, 'Faces', obj.elems(:,obj.sf.contour),'FaceColor','none','EdgeColor','k');
            elseif nargin == 4
                dg     = norm( max(nodes) - min(nodes) );
                maxs = max( abs(min(min(varargin{1}))), abs(max(max(varargin{1})) ) );
                defnodes = (nodes + varargin{1} ./ maxs * dg * varargin{2});
                patch('Vertices', defnodes, 'Faces', obj.elems(:,obj.sf.contour),'FaceColor','none','EdgeColor','r');
            end
        end
        function plot(obj,nodes)
            hold on, axis off;
            daspect([1 1 1]);
            patch('Vertices', nodes, 'Faces', obj.elems(:,obj.sf.contour),'FaceColor','none','EdgeColor','k');
            patch('Vertices', nodes, 'Faces', obj.elems(:,obj.sf.contour),'FaceColor',[0.8 0.8 0.8]);
        end
        function plotMap(obj,nodes,q,C,scd)
            hold, axis off;
            daspect([1 1 1]);
            colormap('jet');
            colorbar;
            nVertices = size(nodes,1);
            if isscalar(C)
                hasNodalResults = isfield(obj.results,'nodal') && isfield(obj.results.nodal,'all') ...
                    && ~isempty(obj.results.nodal.all);
                if hasNodalResults && C == round(C) && C >= 1 && C <= size(obj.results.nodal.all,2)
                    C = obj.results.nodal.all(:,C);
                else
                    error("plotMap expects one color per vertex or a valid nodal result index.");
                end
            elseif isvector(C)
                C = C(:);
            elseif size(C,1) == 3 && size(C,2) == nVertices
                C = C';
            end

            if size(C,1) ~= nVertices
                error("plotMap expected %d vertex colors but received %d.", nVertices, size(C,1));
            end

            patch('Vertices', nodes+scd*q, 'Faces', obj.elems(:,obj.sf.contour), 'FaceVertexCData', C , "FaceColor", "interp", "EdgeColor","none", "FaceAlpha", 1 );
        end
    end
end
