classdef SolidElasticElem < FiniteElement
    methods
        function obj = SolidElasticElem(sf,p)
             obj = obj@FiniteElement(sf,p);
             obj.ndofs=["ux" "uy" "uz"];
             obj.results.names  = ["exx" "eyy" "ezz" "exy" "eyz" "exz" "sxx" "syy" "szz" "sxy" "syz" "sxz" "sHM" "rho"];
             obj.results.descriptions  = ["strain member exx" "strain member eyy" "strain member ezz" ...
                 "strain member exy" "strain member eyz" "strain member exz" "stress member sxx" "stress member syy" "stress member szz" ...
                 "stress member sxy" "stress member syz" "stress member sxz" "Huber-Mises stress" "Top opt density"];
        end
        function setIsotropicMaterial( obj, E, nu, rho )
            D = E / ( 1.0 + nu ) / ( 1 - 2.0 * nu ) * [ 1-nu nu nu 0 0 0; ...
                nu 1-nu nu 0 0 0; nu nu 1-nu 0 0 0;  0 0 0 (1.0 - 2.0 * nu) / 2.0 0 0; ...
                0 0 0 0 (1.0 - 2.0 * nu) / 2.0 0; 0 0 0 0 0 (1.0 - 2.0 * nu) / 2.0 ];
            M =  diag([rho rho rho]); 
            obj.props.D = D;
            obj.props.M = M;
        end
        function B = strainDerivMatrix( ~, dNx )
            ndofs = 3 * size(dNx,3); 
            nip = size(dNx,3); 
            nnd = size(dNx,2); 
            B = cell(nip,1);
            Bn = zeros(6,ndofs);
            for i=1:nip
                for j = 1:nnd
                  Bn(1, 3*j-2) = dNx(1,j,i);
                  Bn(2, 3*j-1) = dNx(2,j,i);
                  Bn(3, 3*j)   = dNx(3,j,i);

                  Bn(4, 3*j-2) = dNx(2,j,i);
                  Bn(5, 3*j-1) = dNx(3,j,i);
                  Bn(6, 3*j-2) = dNx(3,j,i);

                  Bn(4, 3*j-1) = dNx(1,j,i);
                  Bn(5, 3*j)   = dNx(2,j,i);
                  Bn(6, 3*j)   = dNx(1,j,i);
                end
                B{i}=Bn;
            end
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
            %dNx = zeros(size(dN,2),size(dN,1), nip );
            K = zeros( dim , dim, nelems );
            B = zeros(6,dim);
            D = obj.mat.D;
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                Ke = zeros( dim , dim );
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
                    Ke = Ke + abs(detJ) * integrator.weights(i) * B'*D*B;
                end
                K(:,:,k) = x(k)*Ke;
            end
            K=K(:);
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
%                temp=obj.props.nodalTemp(obj.elems(Telems(k),:));
             %   tempCoeff=temp*N;
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
                    %Q=[alpha*tempCoeff(i),alpha*tempCoeff(i),alpha*tempCoeff(i),0,0,0]';
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
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            dNx = zeros(size(dN,2),size(dN,1), nip );
            K = zeros( dim , dim, nelems );
            B = zeros(6,dim);
            dK = zeros( dim, nelems, nsens );
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                for s=1:nsens
                    Ke = zeros( dim , dim );
                    dD = obj.mat.dD(:,:,s);
                    for i=1:nip
                        J = dNtrc{i}*elemX;
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
                        Ke = Ke + abs(detJ) * integrator.weights(i) * B'*dD*B;
                    end
                    dK(:,k,s) = Ke*qelems(:,k);
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
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNx = zeros(size(dN,2),size(dN,1), nip );
            K = zeros( dim , dim, nelems );
            B = zeros(6,dim);
            elemX = nodes(obj.elems(1,:),:);
            Ke = zeros( dim , dim );
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
                Ke = Ke + abs(detJ) * integrator.weights(i) * B'*obj.mat.D*B;
            end
            for k=1:nelems
                K(:,:,k)=x(k)*Ke;
            end
            K=K(:);
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
            obj.results.gp.all(14,:,:) = repmat(x,1,nip);
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
                patch('Vertices', nodes, 'Faces', obj.elems(:,obj.sf.contour),'FaceColor','none','EdgeColor','k');
            elseif nargin == 4
                dg     = norm( max(nodes) - min(nodes) );
                maxs = max( abs(min(min(varargin{1}))), abs(max(max(varargin{1})) ) );
                defnodes = (nodes + varargin{1} ./ maxs * dg * varargin{2});
                allfaces = reshape(obj.elems(:,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems,1))';
                [~,ifaces] = unique( sort(allfaces,2),'rows' );
                patch('Vertices', defnodes, 'Faces', allfaces(ifaces,:),'FaceColor','none','EdgeColor','r');
            end
        end
        function plot(obj,nodes,varargin)
            hold on, axis off;
            daspect([1 1 1]);
            if nargin == 2
                col=[0.8 0.8 0.8];
            else
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
            patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',col,'EdgeColor','k');
            %patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',col);
            %patch('Vertices', nodes, 'Faces', plotfaces,'FaceColor',col,"FaceAlpha",0.3);
        end
        function plotSolidSelected(obj,nodes,elem_inds,varargin)
            hold on, axis off;
            daspect([1 1 1]);
            if nargin == 3
                col=[0.8 0.8 0.8];
            else
                col=varargin{1};
            end
            allfaces = reshape(obj.elems(elem_inds,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(find(elem_inds),1))';
            [~,ifaces] = unique( sort(reshape(obj.elems(elem_inds,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(find(elem_inds),1))',2),'rows' );
            patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor','none','EdgeColor','k');
            patch('Vertices', nodes, 'Faces', allfaces(ifaces,:),'FaceColor',col);
        end
        function plotMap(obj,nodes,q,valueName,scd)
            hold, axis off;
            daspect([1 1 1]);
            colormap('jet');
            colorbar;
            valueIndex = find(obj.results.names == valueName);
            if size(valueIndex,2)==0
                valueIndex = find(obj.ndofs == valueName);
                if size(valueIndex,2)==0
                    error("Map name " + valueName + " not implemented in element:" + class(obj));
                else
                    C = q(:,valueIndex);
                    title("displacement "+valueName);
                end
            else
                C = obj.results.nodal.all(:,valueIndex);
                title(obj.results.descriptions(valueIndex));
            end
             allfaces = reshape(obj.elems(:,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems,1))';
            [~,ifaces] = unique( sort(reshape(obj.elems(:,obj.sf.fcontours)',size(obj.sf.fcontours,1),size(obj.sf.fcontours,2)*size(obj.elems,1))',2),'rows' );
            patch('Vertices', nodes+scd*q, 'Faces', allfaces(ifaces,:), 'FaceVertexCData', C , "FaceColor", "interp", "EdgeColor","none", "FaceAlpha", 1 );
        end
    end
end

