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
        function K = computeStifnessMatrix(obj, nodes, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            if ( nargin == 3 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            K = zeros( dim , dim, nelems );
            B = zeros(3,dim);
            weights = integrator.weights;
            D = obj.mat.D;
            h = obj.props.h;
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                Ke = zeros( dim , dim );
                for i=1:nip
                    J = dNtrc{i}*elemX;
                    detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1);
                    dNx = (1 / detJ * [ J(2,2) -J(1,2); -J(2,1)  J(1,1) ]) * dNtrc{i};
                    for j = 1:nnd
                      B(1, 2*j-1) = dNx(1,j);
                      B(2, 2*j)   = dNx(2,j);
                      B(3, 2*j-1) = dNx(2,j);
                      B(3, 2*j)   = dNx(1,j);
                    end
                    Ke = Ke + abs(detJ) * weights(i) * h * B'*D*B;
                end
                K(:,:,k) = x(k)*Ke;
            end
            K=K(:);
        end
        function K = computeGeometricStifnessMatrix(obj, nodes, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            if ( nargin == 3 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            K = zeros( dim , dim, nelems );
            G = zeros(2,dim);
            B = zeros(3,dim);
            weights = integrator.weights;
            S = zeros(3,3);
            s = zeros(2,2);
            h = obj.props.h;
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                Ke = zeros( dim , dim );
                for i=1:nip
                    J = dNtrc{i}*elemX;
                    detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1);
                    dNx = (1 / detJ * [ J(2,2) -J(1,2); -J(2,1)  J(1,1) ]) * dNtrc{i};
                    S(1,1)=obj.results.gp.stress(1,k,i);
                    S(2,2)=obj.results.gp.stress(2,k,i);
                    S(3,3)=obj.results.gp.stress(3,k,i);
                    
                    s(1,1)=obj.results.gp.stress(1,k,i);
                    s(2,2)=obj.results.gp.stress(2,k,i);
                    s(2,1)=obj.results.gp.stress(3,k,i);
                    s(1,2)=obj.results.gp.stress(3,k,i);
                    for j = 1:nnd
                      G(1, 2*j-1) = dNx(1,j);
                      G(2, 2*j)   = dNx(2,j);

                      B(1, 2*j-1) = dNx(1,j);
                      B(2, 2*j)   = dNx(2,j);
                      B(3, 2*j-1) = dNx(2,j);
                      B(3, 2*j)   = dNx(1,j);
                    end
                    %Ke = Ke + abs(detJ) * weights(i) * h * B'*S*B;
                    Ke = Ke + abs(detJ) * weights(i) * h * G'*s*G;
                end
                K(:,:,k) = x(k)*Ke;
            end
            K=K(:);
        end
        function dK = computeStifnessMatrixGradMat(obj, nodes, q, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            nsens = size(obj.mat.dD,3);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            if ( nargin == 4 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            dK = zeros( dim, nelems, nsens );
            B = zeros(3,dim);
            weights = integrator.weights;
            h = obj.props.h;
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                for s=1:nsens
                    dD = obj.mat.dD(:,:,s);
                    Ke = zeros( dim , dim );
                    for i=1:nip
                        J = dNtrc{i}*elemX;
                        detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1);
                        dNx = (1 / detJ * [ J(2,2) -J(1,2); -J(2,1)  J(1,1) ]) * dNtrc{i};
                        for j = 1:nnd
                          B(1, 2*j-1) = dNx(1,j);
                          B(2, 2*j)   = dNx(2,j);
                          B(3, 2*j-1) = dNx(2,j);
                          B(3, 2*j)   = dNx(1,j);
                        end
                        Ke = Ke + abs(detJ) * weights(i) * h * B'*dD*B;
                    end
                    dK(:,k,s) = x(k)*Ke*qelems(:,k);
                end
            end
            dK=reshape(dK,[nnodes * ndofs*nelems nsens]);
        end
        function dK = computeStifnessMatrixGradX(obj, nodes, q, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            nsens = size(obj.mat.dD,3);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            if ( nargin == 4 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            dK = zeros( dim, nelems );
            B = zeros(3,dim);
            weights = integrator.weights;
            h = obj.props.h;
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            D=obj.mat.D;
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                Ke = zeros( dim , dim );
                for i=1:nip
                    J = dNtrc{i}*elemX;
                    detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1);
                    dNx = (1 / detJ * [ J(2,2) -J(1,2); -J(2,1)  J(1,1) ]) * dNtrc{i};
                    for j = 1:nnd
                      B(1, 2*j-1) = dNx(1,j);
                      B(2, 2*j)   = dNx(2,j);
                      B(3, 2*j-1) = dNx(2,j);
                      B(3, 2*j)   = dNx(1,j);
                    end
                    Ke = Ke + abs(detJ) * weights(i) * h * B'*D*B;
                end
                dK(:,k) = Ke*qelems(:,k);
            end
        end
        function K = computeStifnessMatrixConst(obj, nodes, x)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            K = zeros( dim , dim, nelems );
            B = zeros(3,dim);
            weights = integrator.weights;
            D = obj.mat.D;
            h = obj.props.h;
            elemX = nodes(obj.elems(1,:),:);
            Ke = zeros( dim , dim );
            for i=1:nip
                J = dNtrc{i}*elemX;
                detJ = J(1,1)*J(2,2)-J(1,2)*J(2,1);
                dNx = (1/detJ*[ J(2,2) -J(1,2); -J(2,1)  J(1,1) ])*dNtrc{i};
                for j = 1:nnd
                  B(1, 2*j-1) = dNx(1,j);
                  B(2, 2*j)   = dNx(2,j);
                  B(3, 2*j-1) = dNx(2,j);
                  B(3, 2*j)   = dNx(1,j);
                end
                Ke=Ke+abs(detJ)*weights(i)*h*B'*D*B;
            end
            for k=1:nelems    
                K(:,:,k) = x(k)*Ke;
            end
            K=K(:);
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
            patch('Vertices', nodes+scd*q, 'Faces', obj.elems(:,obj.sf.contour), 'FaceVertexCData', C , "FaceColor", "interp", "EdgeColor","none", "FaceAlpha", 1 );
        end
    end
end

