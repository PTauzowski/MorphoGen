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
            N(1,1:2:nd*nnodes-1,:) = Nsf;
            N(2,2:2:nd*nnodes,:) = Nsf;
        end
        function [J, J1, detJ] = computeJacobian(obj,nodes,dN,el_idx)
            nelems = size(obj.elems,1);
            if (~isempty(el_idx))
                nelems=numel(el_idx);
                elem_nodes = nodes(obj.elems(el_idx,:)',:);
            else
                elem_nodes = nodes(obj.elems',:);
            end
            nnodes = size(obj.elems,2);

            %elemX=permute(repmat(reshape(elem_nodes,nnodes,nelems,2),[1,1,1,4]),[1,3,2,4]); 

            elemX = reshape(elem_nodes, nnodes, nelems, 2);
            elemX = repmat( elemX, [1, 1, 1, 4]);  % then align for pagemtimes
            elemX = permute(elemX, [1, 3, 2, 4]); % final shape: nnodes × 2 × nelems × 4

            J=pagemtimes(dN,elemX); 
            detJ = J(1,1,:,:) .* J(2,2,:,:) - J(1,2,:,:) .* J(2,1,:,:);
            J1 = 1 ./ detJ .* [ J(2,2,:,:) -J(1,2,:,:); -J(2,1,:,:)  J(1,1,:,:) ];
        end
        function B = computeStrainDerivativesMatrix(obj,dNx,nip,el_idx)
            nelems = size(obj.elems,1);
            if (~isempty(el_idx))
                nelems=numel(el_idx);
            end
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            B = zeros(3,dim,nelems,nip);
            cols = 1:2:(2*nnodes);      
            cols2 = cols + 1;       
            
            B(1, cols,   :, :) = dNx(1, :, :, :);
            B(2, cols2,  :, :) = dNx(2, :, :, :);
            B(3, cols,   :, :) = dNx(2, :, :, :);
            B(3, cols2,  :, :) = dNx(1, :, :, :);
        end

        function s=computeGeometricStressMatrix(obj,nip)
            s = zeros(2,2,size(obj.elems,1),nip);
            s(1,1,:,:)=obj.results.gp.stress(1,:,:);
            s(2,2,:,:)=obj.results.gp.stress(2,:,:);
            s(2,1,:,:)=obj.results.gp.stress(3,:,:);
            s(1,2,:,:)=obj.results.gp.stress(3,:,:);
        end

        function K=composeGeometricStifnessMatrix(obj,K,So)
             dim=size(K,1);
             K(1:2:dim,1:2:dim,:) = So;
             K(2:2:dim,2:2:dim,:) = So;
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

        function plotWithSettings(obj, nodes, varargin)
            hold on;
            daspect([1 1 1]);
            nodesplot=nodes;
            plotWired=false;
            plotNodes=false;
            elem_color=[0.8 0.8 0.8];
            edge_color='k';
            nodes_color='r';
            elem_num=1:size(obj.elems,1);
            for k=1:nargin-2
                if isstring(varargin{k})
                    if varargin{k}=="deformed"
                        dg     = norm( max(nodes) - min(nodes) );
                        maxs = max( abs(min(min(varargin{k+1}))), abs(max(max(varargin{k+1})) ) );
                        nodesplot = (nodes + varargin{k+1} ./ maxs * dg * varargin{k+2});
                    end
                    if varargin{k}=="map"
                        C=varargin{k+1};
                    end
                    if varargin{k}=="edge color"
                        edge_color=varargin{k+1};
                    end
                    if varargin{k}=="elem color"
                        elem_color=varargin{k+1};
                    end
                    if varargin{k}=="elem nums"
                        elem_num=varargin{k+1};
                    end
                    if varargin{k}=="wired"
                        plotWired=varargin{k+1};
                        elem_color='r';
                    end
                    if varargin{k}=="nodes"
                        plotNodes=varargin{k+1};
                        nodes_color='r';
                    end
                end
            end
            if exist('C','var')
                colormap('jet');
                colorbar;
                patch('Vertices', nodesplot, 'Faces', obj.elems(elem_num,obj.sf.contour), 'FaceVertexCData', C , "FaceColor", "interp", "EdgeColor","none", "FaceAlpha", 1 );
            else 
                patch('Vertices', nodesplot, 'Faces', obj.elems(elem_num,obj.sf.contour),'FaceColor',elem_color,'EdgeColor',edge_color);
            end
            if plotNodes
                scatter(nodesplot(:,1),nodesplot(:,2),".")
            end
        end
    end
end

