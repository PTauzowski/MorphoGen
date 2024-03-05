classdef PlaneStressElem < PlaneElem
    
   methods
        function obj = PlaneStressElem(sf,elems)
             obj = obj@PlaneElem(sf,elems);
             obj.ndofs=["ux" "uy"];
             obj.results.names  = ["exx" "eyy" "exy" "sxx" "syy" "sxy" "e1" "e2" "gmax" "etheta" "tr(e)" "Vol" "s1" "s2" "tmax" "stheta" "sHM" "rho"];
             obj.results.descriptions  = ["strain member exx" "strain member eyy" "strain member exy" "stress member sxx"...
                 "stress member syy" "stress member sxy" "principal strain e1" "principal strain e2" "maximal shear strain"...
                 "principal strain angle" "strain tensor trace Tr(e)" "volume change" "principal stress s1" "principal stress s2" ...
                 "maximal shear stress" "principal stress angle" "Huber-Mises stress" "Top opt density"];
            
        end
        function setIsotropicMaterial( obj, E, nu, rho )
            D = E/(1-nu*nu)*[ 1  nu 0; ...
                              nu 1  0; ...
                              0  0  (1-nu)/2 ];
            M =  diag([rho rho]); 
            obj.props.D = D;
            obj.props.M = M;
        end
        function initializeResults(obj)
            nelems = size(obj.elems,1);
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            obj.results.gp.strain = zeros(3,nelems,nip);
            obj.results.gp.stress = zeros(3,nelems,nip);
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
            B = zeros(3,dim);
            D = obj.mat.D;
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
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
                    e = B*qelems(:,k); %reshape( q(obj.elems(k,:),:)', nnodes * ndofs,1 );
                    s = x(k)*D*e;
                    obj.results.gp.strain(:,k,i) = e;
                    obj.results.gp.stress(:,k,i) = s;
                end
            end
            %obj.results.gp.strain = permute(strain,[2,3,1]);
            %obj.results.gp.stress = permute(stress,[2,3,1]);
            exx = obj.results.gp.strain(1,:,:);
            eyy = obj.results.gp.strain(2,:,:);
            exy = obj.results.gp.strain(3,:,:);
            sxx = obj.results.gp.stress(1,:,:);
            syy = obj.results.gp.stress(2,:,:);
            sxy = obj.results.gp.stress(3,:,:);
            e1 =  ( exx + eyy ) ./ 2.0 + sqrt( ( (exx - eyy) ./ 2.0 ) .* ( (exx - eyy) ./ 2.0 ) + ( exy .* exy )  );
            e2 =  ( exx + eyy ) ./ 2.0 - sqrt( ( (exx - eyy) ./ 2.0 ) .* ( (exx - eyy) ./ 2.0 ) + ( exy .* exy )  );
            maxt = ( e1 - e2 ) ./ 2.0;
            etheta = 2 .* atan( exy ./ (exx - eyy));
            etr  = exx + eyy;
            vol  = ( 1.0 + e1 ) .* ( 1.0 + e2 );
            s1 =  ( sxx + syy ) ./ 2.0 + sqrt( ( (sxx - syy) ./ 2.0 ) .* ( (sxx - syy) ./ 2.0 ) + ( sxy .* sxy )  );
            s2 =  ( sxx + syy ) ./ 2.0 - sqrt( ( (sxx - syy) ./ 2.0 ) .* ( (sxx - syy) ./ 2.0 ) + ( sxy .* sxy )  );
            maxs = ( s1 - s2 ) ./ 2.0;
            stheta = 2 .* atan( sxy ./ (sxx - syy) );
            sHM = sqrt(  s1 .* s1 - s1 .* s2 + s2 .* s2 );
           
            obj.results.gp.all(1,:,:) = exx(1,:,:);
            obj.results.gp.all(2,:,:) = eyy(1,:,:);
            obj.results.gp.all(3,:,:) = exy(1,:,:);
            obj.results.gp.all(4,:,:) = sxx(1,:,:);
            obj.results.gp.all(5,:,:) = syy(1,:,:);
            obj.results.gp.all(6,:,:) = sxy(1,:,:);
            obj.results.gp.all(7,:,:) = e1(1,:,:); 
            obj.results.gp.all(8,:,:) = e2(1,:,:);
            obj.results.gp.all(9,:,:) = maxt(1,:,:);
            obj.results.gp.all(10,:,:) = etheta(1,:,:);
            obj.results.gp.all(11,:,:) = etr(1,:,:);
            obj.results.gp.all(12,:,:) = vol(1,:,:);
            obj.results.gp.all(13,:,:) = s1(1,:,:);
            obj.results.gp.all(14,:,:) = s2(1,:,:);
            obj.results.gp.all(15,:,:) = maxs(1,:,:);
            obj.results.gp.all(16,:,:) = stheta(1,:,:);
            obj.results.gp.all(17,:,:) = sHM(1,:,:);
            obj.results.gp.all(18,:,:) = repmat(x,1,nip);
        end
        function [HMs, dHMs] = computeHMstress(obj,nodes, nelem, q, ddq, penalty, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            nsens = size(ddq,2);
            dim = nnodes * ndofs;
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            if ( nargin == 7 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            B = zeros(3,dim);
            stress = zeros(3,nsens,nip);
            dstress = zeros(3,nsens,nip);
            %obj.results.GPvalues = zeros(size(obj.results.names,2),nelems,nip);
            D = obj.mat.D;
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            for k=1:nsens
                elemX = nodes(obj.elems(nelem,:),:);
                dqe=reshape(ddq(:,k),ndofs,size(nodes,1))';
                dqelems = reshape( dqe( obj.elems',:)', nnodes * ndofs, nelems );
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
                    e = B*qelems(:,nelem); 
                    s = x(nelem)^penalty*D*e;
                    de = B*dqelems(:,nelem); 
                    if nelem == k
                        ds = x(nelem).^penalty*D*de+penalty*x(nelem).^(penalty-1)*D*e;
                    else
                        ds = x(nelem).^penalty*D*de;
                    end
                    stress(:,k,i) = s;
                    dstress(:,k,i) = ds;
                end
            end
            stress = permute(stress,[2,3,1]);
            dstress = permute(dstress,[2,3,1]);
           
            sxx = stress(:,:,1);
            syy = stress(:,:,2);
            sxy = stress(:,:,3);

            dsxx = dstress(:,:,1);
            dsyy = dstress(:,:,2);
            dsxy = dstress(:,:,3);

            ssqr = sqrt( ( (sxx - syy) ./ 2.0 ) .* ( (sxx - syy) ./ 2.0 ) + ( sxy .* sxy )  );
            dssqr = (sxx - syy).*(dsxx - dsyy)/2 + 2 * sxy .* dsxy;
          
            s1 =  ( sxx + syy ) ./ 2.0 + ssqr;
            s2 =  ( sxx + syy ) ./ 2.0 - ssqr;

            ds1 =  ( dsxx + dsyy ) ./ 2.0 + dssqr./(2*ssqr);
            ds2 =  ( dsxx + dsyy ) ./ 2.0 - dssqr./(2*ssqr);
           
            HMs = sqrt(  s1 .* s1 - s1 .* s2 + s2 .* s2 );
            dHMs = (2*s1 .* ds1 - (ds1 .* s2 + s1 .* ds2) + 2*s2 .* ds2 )./(2*HMs);
           
        end
        
    end
end

