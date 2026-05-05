classdef SolidElasticElem < FiniteElement
    properties
        face_alpha, face_color, edge_color
    end
    methods
        function obj = SolidElasticElem(sf,p)
             obj = obj@FiniteElement(sf,p);
             obj.ndofs=["ux" "uy" "uz"];
             obj.results.names  = ["exx" "eyy" "ezz" "exy" "eyz" "exz" "sxx" "syy" "szz" "sxy" "syz" "sxz" "sHM" "rho" "T"];
             obj.results.descriptions  = ["strain member exx" "strain member eyy" "strain member ezz" ...
                 "strain member exy" "strain member eyz" "strain member exz" "stress member sxx" "stress member syy" "stress member szz" ...
                 "stress member sxy" "stress member syz" "stress member sxz" "Huber-Mises stress" "Top opt density" "Nodal temperature"];
             obj.face_alpha = 1.0;
             obj.edge_color='k';
             obj.face_color=[0.8 0.8 0.8];
        end
        function setIsotropicMaterial( obj, E, nu, rho )
            D = E / ( 1.0 + nu ) / ( 1 - 2.0 * nu ) * [ 1-nu nu nu 0 0 0; ...
                nu 1-nu nu 0 0 0; nu nu 1-nu 0 0 0;  0 0 0 (1.0 - 2.0 * nu) / 2.0 0 0; ...
                0 0 0 0 (1.0 - 2.0 * nu) / 2.0 0; 0 0 0 0 0 (1.0 - 2.0 * nu) / 2.0 ];
            M =  diag([rho rho rho]); 
            obj.props.D = D;
            obj.props.M = M;
        end
        function [Jinv, detJ] = jacobianInversePages(obj, nodes, elemIds, dNtr)
            elemX = obj.elementNodePages(nodes, elemIds);
            nelems = numel(elemIds);
            nip = size(dNtr, 3);
            detJ = zeros(1, 1, nelems, nip);
            Jinv = zeros(3, 3, nelems, nip);
            for ip = 1:nip
                J = pagemtimes(dNtr(:, :, ip), elemX);
                detJ(:, :, :, ip) = J(1,1,:).*J(2,2,:).*J(3,3,:) ...
                    - J(1,2,:).*J(2,1,:).*J(3,3,:) ...
                    - J(1,1,:).*J(2,3,:).*J(3,2,:) ...
                    + J(1,3,:).*J(2,1,:).*J(3,2,:) ...
                    + J(1,2,:).*J(2,3,:).*J(3,1,:) ...
                    - J(1,3,:).*J(2,2,:).*J(3,1,:);
                dj = detJ(:, :, :, ip);
                Jinv(1,1,:,ip) =  (J(2,2,:).*J(3,3,:) - J(2,3,:).*J(3,2,:)) ./ dj;
                Jinv(1,2,:,ip) = -(J(1,2,:).*J(3,3,:) - J(1,3,:).*J(3,2,:)) ./ dj;
                Jinv(1,3,:,ip) =  (J(1,2,:).*J(2,3,:) - J(1,3,:).*J(2,2,:)) ./ dj;
                Jinv(2,1,:,ip) = -(J(2,1,:).*J(3,3,:) - J(2,3,:).*J(3,1,:)) ./ dj;
                Jinv(2,2,:,ip) =  (J(1,1,:).*J(3,3,:) - J(1,3,:).*J(3,1,:)) ./ dj;
                Jinv(2,3,:,ip) = -(J(1,1,:).*J(2,3,:) - J(1,3,:).*J(2,1,:)) ./ dj;
                Jinv(3,1,:,ip) =  (J(2,1,:).*J(3,2,:) - J(2,2,:).*J(3,1,:)) ./ dj;
                Jinv(3,2,:,ip) = -(J(1,1,:).*J(3,2,:) - J(1,2,:).*J(3,1,:)) ./ dj;
                Jinv(3,3,:,ip) =  (J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)) ./ dj;
            end
        end
        function B = strainBPages(obj, Jinv, dNtr)
            nnodes = size(obj.elems, 2);
            dim = 3 * nnodes;
            nelems = size(Jinv, 3);
            nip = size(Jinv, 4);
            B = zeros(6, dim, nelems, nip);
            for ip = 1:nip
                dNx = pagemtimes(Jinv(:, :, :, ip), dNtr(:, :, ip));
                cols = 1:nnodes;
                c1 = 3 * cols - 2;
                c2 = 3 * cols - 1;
                c3 = 3 * cols;
                B(1, c1, :, ip) = dNx(1, cols, :);
                B(2, c2, :, ip) = dNx(2, cols, :);
                B(3, c3, :, ip) = dNx(3, cols, :);
                B(4, c2, :, ip) = dNx(3, cols, :);
                B(4, c3, :, ip) = dNx(2, cols, :);
                B(5, c1, :, ip) = dNx(3, cols, :);
                B(5, c3, :, ip) = dNx(1, cols, :);
                B(6, c1, :, ip) = dNx(2, cols, :);
                B(6, c2, :, ip) = dNx(1, cols, :);
            end
        end
        function Sg = geometricGradientPages(~, Jinv, dNtr)
            nelems = size(Jinv, 3);
            nip = size(Jinv, 4);
            nnodes = size(dNtr, 2);
            Sg = zeros(3, nnodes, nelems, nip);
            for ip = 1:nip
                Sg(:, :, :, ip) = pagemtimes(Jinv(:, :, :, ip), dNtr(:, :, ip));
            end
        end
        function s = geometricStressPages(obj, elemIds, nip)
            nelems = numel(elemIds);
            s = zeros(3, 3, nelems, nip);
            stress = obj.results.gp.stress(elemIds, :, :);
            s(1,1,:,:) = reshape(stress(:,:,1), 1, 1, nelems, nip);
            s(2,2,:,:) = reshape(stress(:,:,2), 1, 1, nelems, nip);
            s(3,3,:,:) = reshape(stress(:,:,3), 1, 1, nelems, nip);
            s(2,3,:,:) = reshape(stress(:,:,4), 1, 1, nelems, nip);
            s(3,2,:,:) = reshape(stress(:,:,4), 1, 1, nelems, nip);
            s(1,3,:,:) = reshape(stress(:,:,5), 1, 1, nelems, nip);
            s(3,1,:,:) = reshape(stress(:,:,5), 1, 1, nelems, nip);
            s(1,2,:,:) = reshape(stress(:,:,6), 1, 1, nelems, nip);
            s(2,1,:,:) = reshape(stress(:,:,6), 1, 1, nelems, nip);
        end
        function K = composeGeometricPages(~, So, scale)
            nnodes = size(So, 1);
            nelems = size(So, 3);
            dim = 3 * nnodes;
            K = zeros(dim, dim, nelems);
            So = So .* reshape(scale, 1, 1, []);
            K(1:3:dim, 1:3:dim, :) = So;
            K(2:3:dim, 2:3:dim, :) = So;
            K(3:3:dim, 3:3:dim, :) = So;
        end
        function N = shapeMatrix( obj, points )
            nnodes = size(obj.elems, 2 );
            nd = size( obj.ndofs, 2 );
            np = size( points, 1 );
            N = zeros( nd, nd*nnodes, np );
            Nsf = obj.sf.computeValue(points);
            for k=1:np
                N(1,1:3:nd*nnodes-2,k) = Nsf(k,:);
                N(2,2:3:nd*nnodes-1,k) = Nsf(k,:);
                N(3,3:3:nd*nnodes,k)   = Nsf(k,:);
            end
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
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [Jinv, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                B = obj.strainBPages(Jinv, dNtr);
                K(:, :, elemIds) = obj.integratePagematrix(B, D, detJ, integrator.weights, x(elemIds));
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
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [Jinv, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                Sg = obj.geometricGradientPages(Jinv, dNtr);
                s = obj.geometricStressPages(elemIds, size(integrator.points, 1));
                So = obj.integratePagematrix(Sg, s, detJ, integrator.weights, ones(numel(elemIds), 1));
                K(:, :, elemIds) = obj.composeGeometricPages(So, x(elemIds));
            end
            K = obj.flattenElementMatrices(K);
        end
        function M = computeMassMatrix(obj, nodes, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            N = obj.shapeMatrix(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            x = obj.elementScale(nelems, varargin{:});
            M = zeros(dim, dim, nelems);
            chunkSize = obj.assemblyChunkSize(dim);
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
                    integrator.weights, x(elemIds));
            end
            M = obj.flattenElementMatrices(M);
        end
        function Pnodal = thermalLoad(obj, nodes, Telems, Pnodal, alpha, varargin)
            nelems = size(Telems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            N=obj.sf.computeValue( integrator.points );
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            if ( nargin == 6 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            %dNx = zeros(size(dN,2),size(dN,1), nip );
            B = zeros(6,dim);
            Q=[alpha,alpha,alpha,0,0,0]';
            obj.props.thermal(Telems)=alpha;
            for k=1:nelems
                elemX = nodes(obj.elems(Telems(k),:),:);
                temps=obj.props.ndT(obj.elems(Telems(k),:));
                tempGP=temps*N';
                Pe = zeros( dim , 1 );
                for i=1:nip
                    J = dNtrc{i}*elemX;
                    detJ = J(1,1)*J(2,2)*J(3,3)-J(1,2)*J(2,1)*J(3,3)-J(1,1)*J(2,3)*J(3,2)+J(1,3)*J(2,1)*J(3,2)+J(1,2)*J(2,3)*J(3,1)-J(1,3)*J(2,2)*J(3,1);
                    invJ   = [ (J(2,2)*J(3,3)-J(2,3)*J(3,2))	-(J(1,2)*J(3,3)-J(1,3)*J(3,2))  (J(1,2)*J(2,3)-J(1,3)*J(2,2) ); ...
              		          -(J(2,1)*J(3,3)-J(2,3)*J(3,1))	 (J(1,1)*J(3,3)-J(1,3)*J(3,1)) -(J(1,1)*J(2,3)-J(1,3)*J(2,1) ); ...
              		           (J(2,1)*J(3,2)-J(2,2)*J(3,1))	-(J(1,1)*J(3,2)-J(1,2)*J(3,1))  (J(1,1)*J(2,2)-J(1,2)*J(2,1) ) ]/detJ;
                    dNx = invJ * dNtr(:,:,i);
                    for j = 1:nnd
                          B(1, 3*j-2) = dNx(1,j);
                          B(2, 3*j-1) = dNx(2,j);
                          B(3, 3*j)   = dNx(3,j);

                          B(4, 3*j-1) = dNx(3,j);
                          B(4, 3*j) = dNx(2,j);
                          
                          B(5, 3*j-2) = dNx(3,j);
                          B(5, 3*j) = dNx(1,j);
                          
                          B(6, 3*j-2) = dNx(2,j);
                          B(6, 3*j-1) = dNx(1,j);

                    end  
                    Q=[alpha*tempGP(i),alpha*tempGP(i),alpha*tempGP(i),0,0,0]';
                    Pe = Pe + abs(detJ) * integrator.weights(i) * B'*Q;
                end
                Pnodal(obj.elems(Telems(k),:),:) = Pnodal(obj.elems(Telems(k),:),:) + reshape(Pe,3,27)';
            end
        end
        function dK = computeStifnessMatrixGradMat(obj, nodes, q)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            nsens = size(obj.mat.dD,3);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            dK = zeros(dim, nelems, nsens);
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            chunkSize = obj.assemblyChunkSize(dim);
            for first = 1:chunkSize:nelems
                elemIds = first:min(first + chunkSize - 1, nelems);
                [Jinv, detJ] = obj.jacobianInversePages(nodes, elemIds, dNtr);
                B = obj.strainBPages(Jinv, dNtr);
                qpages = reshape(qelems(:, elemIds), dim, 1, []);
                for s=1:nsens
                    Ke = obj.integratePagematrix(B, obj.mat.dD(:,:,s), detJ, ...
                        integrator.weights, ones(numel(elemIds), 1));
                    dK(:, elemIds, s) = squeeze(pagemtimes(Ke, qpages));
                end
            end
            dK=reshape(dK,[nnodes*ndofs*nelems nsens]);
        end
        function K = computeStifnessMatrixConst(obj, nodes, x)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            dN = obj.sf.computeGradient(integrator.points);
            dNtr = permute(dN,[2,1,3]);
            x = obj.elementScale(nelems, x);
            [Jinv, detJ] = obj.jacobianInversePages(nodes, 1, dNtr);
            B = obj.strainBPages(Jinv, dNtr);
            Ke = obj.integratePagematrix(B, obj.mat.D, detJ, integrator.weights, 1);
            K = Ke .* reshape(x, 1, 1, []);
            K = obj.flattenElementMatrices(K);
        end
        function initializeResults(obj)
            nelems = size(obj.elems,1);
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            obj.results.gp.strain = zeros(6,nelems,nip);
            obj.results.gp.stress = zeros(6,nelems,nip);
            obj.results.gp.all = zeros(size(obj.results.names,2),nelems,nip);
        end
        function computeResults(obj,nodes, q, varargin)
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
            if ( nargin == 4 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            B = zeros(6,dim);
            strain = zeros(6,nelems,nip);
            stress = zeros(6,nelems,nip);
            obj.results.gp.all = zeros(size(obj.results.names,2),nelems,nip);
            D = obj.mat.D;
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            for k=1:nelems
                et=[obj.props.thermal(k), obj.props.thermal(k), obj.props.thermal(k),0,0,0]';
                elemX = nodes(obj.elems(k,:),:);
                for i=1:nip
                    J = dNtr(:,:,i) * elemX;
                    detJ = J(1,1)*J(2,2)*J(3,3)-J(1,2)*J(2,1)*J(3,3)-J(1,1)*J(2,3)*J(3,2)+J(1,3)*J(2,1)*J(3,2)+J(1,2)*J(2,3)*J(3,1)-J(1,3)*J(2,2)*J(3,1);
                    invJ   = [ (J(2,2)*J(3,3)-J(2,3)*J(3,2))/detJ	-(J(1,2)*J(3,3)-J(1,3)*J(3,2))/detJ  (J(1,2)*J(2,3)-J(1,3)*J(2,2) )/detJ; ...
              		          -(J(2,1)*J(3,3)-J(2,3)*J(3,1))/detJ	 (J(1,1)*J(3,3)-J(1,3)*J(3,1))/detJ -(J(1,1)*J(2,3)-J(1,3)*J(2,1) )/detJ; ...
              		           (J(2,1)*J(3,2)-J(2,2)*J(3,1))/detJ	-(J(1,1)*J(3,2)-J(1,2)*J(3,1))/detJ  (J(1,1)*J(2,2)-J(1,2)*J(2,1) )/detJ ];
                    dNx = invJ * dNtr(:,:,i);
                    for j = 1:nnd
                          B(1, 3*j-2) = dNx(1,j);
                          B(2, 3*j-1) = dNx(2,j);
                          B(3, 3*j)   = dNx(3,j);

                          B(4, 3*j-1) = dNx(3,j);
                          B(4, 3*j) = dNx(2,j);
                          
                          B(5, 3*j-2) = dNx(3,j);
                          B(5, 3*j) = dNx(1,j);
                          
                          B(6, 3*j-2) = dNx(2,j);
                          B(6, 3*j-1) = dNx(1,j);
                    end
                    e = B*qelems(:,k); 
                    stress(:,k,i) = x(k)*D*(e-et);
                    strain(:,k,i) = e;
                end
            end
            obj.results.gp.strain = permute(strain,[2,3,1]);
            obj.results.gp.stress = permute(stress,[2,3,1]);
            exx = obj.results.gp.strain(:,:,1);
            eyy = obj.results.gp.strain(:,:,2);
            ezz = obj.results.gp.strain(:,:,3);
            exy = obj.results.gp.strain(:,:,4);
            eyz = obj.results.gp.strain(:,:,5);
            exz = obj.results.gp.strain(:,:,6);
            sxx = obj.results.gp.stress(:,:,1);
            syy = obj.results.gp.stress(:,:,2);
            szz = obj.results.gp.stress(:,:,3);
            sxy = obj.results.gp.stress(:,:,4);
            syz = obj.results.gp.stress(:,:,5);
            sxz = obj.results.gp.stress(:,:,6);
            sHM = sqrt( 1/2*( (sxx-syy).^2+(syy-szz).^2+(szz-sxx).^2 )+3*( sxy.^2+syz.^2+sxz.^2) ) ;
            obj.results.gp.all(1,:,:) = exx;
            obj.results.gp.all(2,:,:) = eyy;
            obj.results.gp.all(3,:,:) = ezz;
            obj.results.gp.all(4,:,:) = exy;
            obj.results.gp.all(5,:,:) = eyz;
            obj.results.gp.all(6,:,:) = exz;
            obj.results.gp.all(7,:,:) = sxx;
            obj.results.gp.all(8,:,:) = syy;
            obj.results.gp.all(9,:,:) = szz;
            obj.results.gp.all(10,:,:) = sxy;
            obj.results.gp.all(11,:,:) = syz;
            obj.results.gp.all(12,:,:) = sxz;
            obj.results.gp.all(13,:,:) = sHM;
            obj.results.gp.all(14,:,:) = repmat(reshape(x(:)', 1, nelems, 1), 1, 1, nip);
        end
        function faces = findFaces( obj, fnodes )
              allfaces = obj.multiObjectList( obj.sf.faces );
              bfaces = false( size(allfaces,1), 1 );
              dfaces = setdiff( allfaces, fnodes ) ;
              for k=1:size(allfaces,1)
                bfaces(k) = isempty( setdiff( allfaces(k,:), fnodes ) );
              end  
              faces = allfaces( bfaces, : );
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
                lPoints= N*elemX;
                dXY = dN * elemX;
                Pfn = zeros( size(elemX,1), size( obj.ndofs,2) );
                lValues = valueFn( lPoints );
                Pfn(:,inds) = lValues
                if mode =="local"
                    nXYZ = dXY./vecnorm(dXY')';
                    tXYZ = [ nXY(:,2) -nXY(:,1) nXY(:,3) ];
                    Pfn = nXYZ .* Pfn(:,1) + tXYZ .* Pfn(:,2) + tXYZ .* Pfn(:,3);
                end
                Pe = zeros( size(Pfn) );
                volume = 0;
                for i=1:np
                    dt = sum( dXY(i,:) .* dXY(i,:) );
                    Pe = Pe + integrator.weights( i ) * sqrt( dt ) * N(i,:)' .* Pfn(i,:);
                    volume = volume + integrator.weights( i ) * sqrt( dt );
                end
                P(edges(k,:),di) = P(edges(k,:),di) + Pe(:,inds);
            end

        end
        function [P, volume] = loadSurfaceIntegral(obj, mode, nodes, faces, dofnames, di, P, valueFn)
            inds = obj.findDofIndices( dofnames );
            integrator = obj.sf.facesf.createIntegrator();
            N = obj.sf.facesf.computeValue( integrator.points );
            dN = obj.sf.facesf.computeGradient( integrator.points );
            nip = size(dN,2);
            np  = size(dN,3);
            for k=1:size(faces,1)
                elemX = nodes(faces(k,:),:);
                Pfn = zeros( size(N,1), size( obj.ndofs,2) );
                Pfn(:,inds) = valueFn( N*elemX );
                if mode =="local"
%                     nXYZ = dXY./vecnorm(dXY')';
%                     tXYZ = [ nXY(:,2) -nXY(:,1) nXY(:,3) ];
%                     Pfn = nXYZ .* Pfn(:,1) + tXYZ .* Pfn(:,2) + tXYZ .* Pfn(:,3);
                end
                Pe = zeros( size(Pfn) );
                volume = 0;
                for i=1:np
                    J = dN(:,:,i)' * elemX;
                    J=J';
                    H =  ( J(1,1) * J(2,2) - J(2,1) * J(1,2) ) * ( J(1,1) * J(2,2) - J(2,1) * J(1,2) ) + ...
                         ( J(2,1) * J(3,2) - J(3,1) * J(2,2) ) * ( J(2,1) * J(3,2) - J(3,1) * J(2,2) ) + ...
                         ( J(3,1) * J(1,2) - J(1,1) * J(3,2) ) * ( J(3,1) * J(1,2) - J(1,1) * J(3,2) );
                    Pe = Pe + integrator.weights( i ) * sqrt( H ) * N(i,:)' .* Pfn(i,:);
                    volume = volume + integrator.weights( i ) * sqrt( H );
                end
                P(faces(k,:),di) = P(faces(k,:),di) + Pe(:,inds);
            end

        end
        function plotWired(obj,nodes,varargin)
            hold on;
            daspect([1 1 1]);
            if nargin==2
                if isempty(obj.selectedElems)
                    patch('Vertices', nodes, 'Faces', obj.elems(:,obj.sf.contour),'FaceColor','none','EdgeColor','k');
                else
                    patch('Vertices', nodes, 'Faces', obj.elems(obj.selectedElems,obj.sf.contour),'FaceColor','none','EdgeColor','k');
                end
            elseif nargin == 4
                dg     = norm( max(nodes) - min(nodes) );
                maxs = max( abs(min(min(varargin{1}))), abs(max(max(varargin{1})) ) );
                defnodes = (nodes + varargin{1} ./ maxs * dg * varargin{2});
                if isempty(obj.selectedElems)
                    allfaces = reshape(obj.elems(:,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems,1))';
                else
                    allfaces = reshape(obj.elems(obj.selectedElems,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems(obj.selectedElems,:),1))';
                end
                [~,ifaces] = unique( sort(allfaces,2),'rows' );
                p=patch('Vertices', defnodes, 'Faces', allfaces(ifaces,:),'FaceColor','none','EdgeColor','k');
                p.LineWidth = 0.01;
            end
        end
        function plot(obj,nodes,varargin)
            col=[0.8 0.8 0.8];
            if nargin > 2
                col=varargin{1};
            end
            allfaces = reshape(obj.elems(:,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems,1))';
            %allfaces = obj.elems(obj.sf.fcontours,:);
            [~,ifaces] = unique( sort(allfaces,2), 'rows' );
            A=allfaces(ifaces,:);
            delfaces=allfaces;
            delfaces(ifaces,:)=[];
            [~,ifaces,~]=setxor( sort(A,2), sort(delfaces,2), 'rows' );
            plotfaces=A(ifaces,:);
            p=patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',obj.face_color,'EdgeColor',obj.edge_color,"FaceAlpha",obj.face_alpha);
            p.LineWidth = 0.01;
            %patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',col);
            %patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',col,"FaceAlpha",0.3);
        end
        function plotSolidDeformed(obj,nodes,qnodal,scale,elem_inds,varargin)
          

            sel_inds=true(size(obj.elems,1),1);
            if ~isempty(obj.selectedElems)
                sel_inds=false(size(obj.elems,1),1);
                sel_inds(obj.selectedElems)=true;
                sel_inds=sel_inds & elem_inds;
            end

            dg  = norm( max(nodes) - min(nodes) );
            maxs = max( abs(min(min(qnodal))), abs(max(max(qnodal)) ) );
            defnodes = (nodes + qnodal ./ maxs * dg * scale);
            allfaces = reshape(obj.elems(sel_inds,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(find(sel_inds),1))';
            %allfaces = obj.elems(obj.sf.fcontours,:);
            [~,ifaces] = unique( sort(allfaces,2), 'rows' );
            A=allfaces(ifaces,:);
            delfaces=allfaces;
            delfaces(ifaces,:)=[];
            [~,ifaces,~]=setxor( sort(A,2), sort(delfaces,2), 'rows' );
            plotfaces=A(ifaces,:);
            patch('Vertices', defnodes, 'Faces', plotfaces,'FaceColor',obj.face_color,'EdgeColor',obj.edge_color,"FaceAlpha",obj.face_alpha)
            %patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',col);
            %patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',col,"FaceAlpha",0.3);
        end
        function plotSolidSelected(obj,nodes,elem_inds,varargin)
            if isempty(elem_inds) || ~any(elem_inds)
                return
            end
            if nargin == 3
                col=[0.8 0.8 0.8];
            else
                col=varargin{1};
            end
            sel_inds=elem_inds;
            if ~isempty(obj.selectedElems)
                sel_inds=false(size(obj.elems,1),1);
                sel_inds(obj.selectedElems)=true;
                sel_inds=sel_inds & elem_inds;
            end
            allfaces = reshape(obj.elems(sel_inds,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(find(sel_inds),1))';
            alledges = reshape(obj.elems(sel_inds,obj.sf.edges)',size(obj.sf.edges,1),size(obj.sf.edges,2)*size(find(sel_inds),1))';
            [~,ifaces,~] = unique( sort(reshape(obj.elems(sel_inds,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(find(sel_inds),1))',2),'rows' );
            [uniqueEdges, ~, idx] = unique(sort(reshape(obj.elems(sel_inds,obj.sf.edges)',size(obj.sf.edges,1),size(obj.sf.edges,2)*size(find(sel_inds),1))',2), 'rows', 'stable');
            counts = histcounts(idx, 1:(max(idx)+1));
            edgesNoDuplicates = uniqueEdges(counts == 1, :);
            patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor',obj.face_color,'EdgeColor',obj.edge_color);
            patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor',col,'EdgeColor','none',"FaceAlpha",obj.face_alpha);
            x=[nodes(edgesNoDuplicates(:,1),1) nodes(edgesNoDuplicates(:,2),1) NaN(size(edgesNoDuplicates,1),1) ];
            y=[nodes(edgesNoDuplicates(:,1),2) nodes(edgesNoDuplicates(:,2),2) NaN(size(edgesNoDuplicates,1),1) ];
            z=[nodes(edgesNoDuplicates(:,1),3) nodes(edgesNoDuplicates(:,2),3) NaN(size(edgesNoDuplicates,1),1) ];
            x = reshape( x', 3 * size(edgesNoDuplicates,1),1);
            y = reshape( y', 3 * size(edgesNoDuplicates,1),1);
            z = reshape( z', 3 * size(edgesNoDuplicates,1),1);
            %line(x,y,z,'Color','k','LineWidth', 2);
            %patch(x,y,z,'EdgeColor','k','Marker','.','MarkerFaceColor','flat');

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

            if isempty(obj.selectedElems)
                allfaces = reshape(obj.elems(:,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems,1))';
                [~,ifaces] = unique( sort(reshape(obj.elems(:,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems,1))',2),'rows' );
                patch('Vertices', nodes+scd*q, 'Faces', allfaces(ifaces,:), 'FaceVertexCData', C , "FaceColor", "interp", "EdgeColor",obj.edge_color, "LineStyle", "none","FaceAlpha",obj.face_alpha);
            else

                allfaces = reshape(obj.elems(obj.selectedElems,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems(obj.selectedElems,:),1))';
                [~,ifaces] = unique( sort(reshape(obj.elems(obj.selectedElems,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems(obj.selectedElems,:),1))',2),'rows' );
                patch('Vertices', nodes+scd*q, 'Faces', allfaces(ifaces,:), 'FaceVertexCData', C , "FaceColor", "interp", "EdgeColor",obj.edge_color, "LineStyle", "none", "FaceAlpha",obj.face_alpha );
            end
        end
        end

end
