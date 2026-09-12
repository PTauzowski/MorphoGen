classdef (Abstract) FiniteElement < handle
   
    properties
       props;
       elems;
       ndofs;
       sf;
       results;
       mat;
       selectedElems;
    end
    
    methods(Abstract)
      N = shapeMatrix( obj, points );
      P = loadLineIntegral(obj, mode, nodes, fedges,dofnames,valueFn);
      initializeResults(obj);
      computeResults(obj,nodes,q,varargin);
      plotWired(nodes,varargin);
      plot(nodes,varargin);
      plotMap(nodes,q,mapName,scd);
    end
    
    methods
        function obj = FiniteElement(sf, elems)
                     obj.sf=sf;
                     obj.elems = elems;
                     obj.props.h=1;
                     obj.props.thermal=zeros(size(elems,1),1);
                     obj.selectedElems=[];
        end
        function setMaterial(obj, mat)
                     obj.mat=mat;
        end
        function i = findDofIndices( obj, sdof ) 
            i = zeros(size(sdof));
            for k=1:size(i,2)
                i(k) = find( obj.ndofs == sdof(k) );
            end
        end
        function multiList = multiObjectList( obj, ElemObjectList )
             multiList = reshape( obj.elems(:,ElemObjectList)',size( ElemObjectList ,1), size( ElemObjectList,2) * size( obj.elems, 1) )';
        end
        function edges = findEdges( obj, fnodes )
            alledges = obj.multiObjectList( obj.sf.edges );
            bedges = false( size(alledges,1), 1 );
            for k=1:size(alledges,1)
                bedges(k) = isempty( setdiff( alledges(k,:), fnodes ) );
            end  
            edges = alledges( bedges, : );
        end
        function  [I,J,V,Ksize] = sparseMatrixAllocDataUniform( obj, gdofs )
            nelems = size( obj.elems, 1 );
            nnodes = size( obj.elems, 2 );
            % Map local element DOFs into the global list without reordering them.
            [tf,idofs] = ismember(obj.ndofs, gdofs);
            assert(all(tf), 'Local DOF missing in global DOF list');
            idofs = reshape(idofs, 1, []);
            ndofs = numel(idofs);
            ngdofs = numel(gdofs);
            Kdim  = ndofs * nnodes;
            Ksize = Kdim * Kdim;
            [ix, iy] = meshgrid( 1:Kdim, 1:Kdim );
            dofOffsets = repmat(idofs, nelems, nnodes);
            alldofs = (repelem(obj.elems, 1, ndofs) - 1) * ngdofs + dofOffsets;
            I = alldofs(1:nelems,ix(:));
            J = alldofs(1:nelems,iy(:));
            V = alldofs;
        end

        function K = computeElementMatrices(obj, scale, weights, detJ, B, D)
            Ke =  reshape( scale , 1, 1, [], 1) .* reshape( weights , 1, 1, 1, []) .* detJ .* pagemtimes(pagemtimes(B,'transpose',D,'none'),B);
            K=sum(Ke,4);
        end

        function P = computeElementSelfWeight(obj, scale, weights, detJ, N, rho)
            Pe =  reshape( scale , 1, [], 1) .* reshape( rho , 1, [], 1) .* reshape( weights , 1, 1, []) .* pagemtimes(detJ , N);
            P=sum(Pe,3);
        end

        function K = computeStifnessMatrix(obj, nodes, el_idx)
            nelems = size(obj.elems,1);
            if (~isempty(el_idx))
                nelems=numel(el_idx);
            end
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = permute(repmat(obj.sf.computeGradient( integrator.points ),[1,1,1,nelems]),[2,1,4,3]);
            [~, J1, detJ] = obj.computeJacobian(nodes,dN,el_idx);
            dNx = pagemtimes(J1,dN);            
            B = obj.computeStrainDerivativesMatrix(dNx,nip,el_idx);
            h = repelem(obj.props.h,nelems,1);
            K = obj.computeElementMatrices(h, integrator.weights, detJ, B, obj.mat.D);
        end

         function M = computeMassMatrix(obj, nodes, el_idx)
            nelems = size(obj.elems,1);
             if (~isempty(el_idx))
                nelems=numel(el_idx);
            end
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = permute(repmat(obj.sf.computeGradient( integrator.points ),[1,1,1,nelems]),[2,1,4,3]);

            [~, ~, detJ] = obj.computeJacobian(nodes,dN,el_idx);
            N = permute(repmat(obj.shapeMatrix( integrator.points ),[1,1,1,nelems]),[1,2,4,3]);
            h = repelem(obj.props.h,nelems,1);
            M = obj.computeElementMatrices(h, integrator.weights, detJ, N, obj.mat.M);
         end


         

         function K = computeGeometricStifnessMatrix(obj, nodes, el_idx)
            nelems = size(obj.elems,1);
            if (~isempty(el_idx))
                nelems=numel(el_idx);
            end
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = permute(repmat(obj.sf.computeGradient( integrator.points ),[1,1,1,nelems]),[2,1,4,3]);

            [~, J1, detJ] = obj.computeJacobian(nodes,dN,el_idx);
            dNx = pagemtimes(J1,dN);            
            
            K = zeros( dim , dim, nelems );
            So = obj.computeElementMatrices( 1, integrator.weights, detJ, dNx, obj.computeGeometricStressMatrix(nip) );
            K = obj.composeGeometricStifnessMatrix(K,So);
         end

         function K = computeGeometricStifnessMatrixOld(obj, nodes, el_idx)
            nelems = size(obj.elems,1);
            if (~isempty(el_idx))
                nelems=numel(el_idx);
            end
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = permute(repmat(obj.sf.computeGradient( integrator.points ),[1,1,1,nelems]),[2,1,4,3]);

            [~, J1, ~] = obj.computeJacobian(nodes,dN,el_idx);
            dNx = pagemtimes(J1,dN);

            s=computeBucklingStressMatrix();

            s = zeros(3,3,nelems,nip);
            s(1,1,:,:)=obj.results.gp.stress(1,:,:);
            s(1,1,:,:)=obj.results.gp.stress(1,:,:);
            s(2,2,:,:)=obj.results.gp.stress(2,:,:);
            s(3,3,:,:)=obj.results.gp.stress(3,:,:);
            s(2,3,:,:)=obj.results.gp.stress(4,:,:);
            s(1,3,:,:)=obj.results.gp.stress(5,:,:);
            s(1,2,:,:)=obj.results.gp.stress(6,:,:);
            s(3,2,:,:)=obj.results.gp.stress(4,:,:);
            s(3,1,:,:)=obj.results.gp.stress(5,:,:);
            s(2,1,:,:)=obj.results.gp.stress(6,:,:);
            K = zeros( dim , dim, nelems );
            So = obj.computeElementMatrices( h, integrator.weights, detJ, dNx, s);
            K=obj.composeGeometricStifnessMatrix(So)
            K(1:3:dim,1:3:dim,k) = x(k)*So;
            K(2:3:dim,2:3:dim,k) = x(k)*So;
            K(3:3:dim,3:3:dim,k) = x(k)*So;

        end
        
        function fromGPToNodal(obj,nnodes)
              GPresults = permute( obj.results.GPvalues,[3,1,2]);
              gpres = GPresults( 1,:,: );  % gp x results
              el = obj.elems;
              nres = zeros( nnodes, size( gpres, 2 ) );
              ires = zeros( nnodes, size( gpres, 2 ) );
              sfv = obj.sf.getRecoveryMatrix();
              for k=1:size(el,1)
                  neres = sfv * GPresults(:,:,k); %tensorprod(sfv,GPresults,3)
                  nres( el( k, : ), : ) = nres( el( k, : ), : ) + neres;
                  ires( el( k, : ), : ) = ires( el( k, : ), : ) + 1;
              end
              obj.results.nodal = nres ./ ires;
        end
        function qelems = createElemSolutionVectors(obj,q)
            qelems = reshape( q( obj.elems',:)', size(obj.elems,2) * size(obj.ndofs,2), size(obj.elems,1) );
        end
        function chunkSize = assemblyChunkSize(obj, matrixDim)
            % Keep paged element matrices below a moderate memory footprint.
            % The default targets roughly 128 MB for the dominant dim x dim
            % page array, leaving room for Jacobians, B-matrices and solver
            % data in parallel workers.
            nelems = size(obj.elems, 1);
            targetBytes = 128 * 1024^2;
            bytesPerElem = max(1, matrixDim * matrixDim) * 8;
            chunkSize = max(1, floor(targetBytes / bytesPerElem));
            chunkSize = min(nelems, max(512, chunkSize));
        end
        function x = elementScale(~, nelems, varargin)
            if isempty(varargin)
                x = ones(nelems, 1);
            else
                x = varargin{1};
                if isscalar(x)
                    x = repmat(x, nelems, 1);
                else
                    x = x(:);
                end
            end
        end
        function K = flattenElementMatrices(~, Kpages)
            % Sparse assemblers in this codebase expect each element matrix
            % flattened in MATLAB column-major page order.
            K = Kpages(:);
        end
        function elemX = elementNodePages(obj, nodes, elemIds)
            nnodes = size(obj.elems, 2);
            spatialDim = size(nodes, 2);
            elemNodes = nodes(obj.elems(elemIds, :)', :);
            elemX = permute(reshape(elemNodes, nnodes, numel(elemIds), spatialDim), [1 3 2]);
        end
        function K = integratePagematrix(~, B, D, detJ, weights, scale)
            % Integrate B' * D * B over all pages. B has size
            % nstrain x ndof x nelem x nip; D can be constant or paged
            % nstrain x nstrain x nelem x nip.
            nelems = size(B, 3);
            dim = size(B, 2);
            K = zeros(dim, dim, nelems);
            for ip = 1:numel(weights)
                Bip = B(:, :, :, ip);
                if ndims(D) <= 2
                    DB = pagemtimes(D, Bip);
                else
                    DB = pagemtimes(D(:, :, :, ip), Bip);
                end
                Kip = pagemtimes(Bip, 'transpose', DB, 'none');
                K = K + abs(detJ(:, :, :, ip)) .* weights(ip) .* Kip;
            end
            K = K .* reshape(scale, 1, 1, []);
        end
        
    end
    
end
