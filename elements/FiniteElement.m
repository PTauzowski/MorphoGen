classdef (Abstract) FiniteElement < handle
   
    properties
       props;
       elems;
       ndofs;
       sf;
       results;
       mat;
    end
    
    methods(Abstract)
      N = shapeMatrix( obj, points );
      P = loadLineIntegral(obj, mode, nodes, fedges, dofnames,valueFn);
      initializeResults(obj);
      computeResults(obj,nodes,q,varargin);
      plotWired(nodes,varargin);
      plot(nodes);
      plotMap(nodes,q,mapName,scd);
    end
    
    methods
        function obj = FiniteElement(sf, elems)
                     obj.sf=sf;
                     obj.elems = elems;
                     obj.props.h=1;
                     obj.props.thermal=zeros(size(elems,1),1);
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
            [~,~,idofs] = intersect(obj.ndofs,gdofs);
            Kdim  = size(obj.ndofs,2) * nnodes;
            Ksize = Kdim * Kdim;
            [ix, iy] = meshgrid( 1:Kdim, 1:Kdim );
            alldofs = (repelem( obj.elems, 1, size(obj.ndofs,2))-1)*size(gdofs,2)+repmat(idofs',nelems,nnodes);
            I = alldofs(1:nelems,ix(:));
            J = alldofs(1:nelems,iy(:));
            V = alldofs;
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
        function K = createShapeBasedElemMatrix(obj, nodes, integrator )
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            nudofs = size( obj.ndofs,2);
            dim = nnodes * nudofs;
            nip = size(integrator.points,1);
            N = obj.shapeMatrix( integrator.points );
            Ntr = permute(N,[2,1,3]);
            dN = obj.sf.computeGradient( integrator.points );
            dK = zeros( dim , dim , nip );
            K = zeros( dim , dim, nelems );
            for k=1:nelems
                [ ~, ~, detJ ] = obj.jacobi( dN, nodes(obj.elems(k,:),:) );
                for i=1:nip
                    dK(:,:,i) = abs(detJ(i)) * integrator.weights(i) * Ntr(:,:,i) * obj.props.M * N(:,:,i);
                end
                K(:,:,k) = sum( dK, 3 );
            end
        end
        function K = createGradBasedElemMatrix(obj, nodes, integrator, matrixB)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            nudofs = size( obj.ndofs,2);
            dim = nnodes * nudofs;
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            dNtr = permute(dN,[2,1,3]);
            dNx = zeros(size(dN,2),size(dN,1), nip );
            K = zeros( dim , dim, nelems );
            for k=1:nelems
                [ ~, invJ, detJ ] = obj.jacobi( dN, nodes(obj.elems(k,:),:) );
                for i=1:nip
                    dNx(:,:,i) = invJ(:,:,i) * dNtr(:,:,i);
                end
                B = matrixB( dNx );
                Ke = zeros( dim , dim );
                for i=1:nip
                    Ke = Ke + abs(detJ(i)) * integrator.weights(i) * B{i}'*obj.props.D*B{i};
                end
                K(:,:,k) = Ke;
            end
        end
    end
    
end

