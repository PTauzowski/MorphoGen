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
        function computeResults(obj,nodes,q, varargin)
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
            strain = zeros(3,nelems,nip);
            stress = zeros(3,nelems,nip);
            obj.results.GPvalues = zeros(size(obj.results.names,2),nelems,nip);
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
                    strain(:,k,i) = e;
                    stress(:,k,i) = s;
                end
            end
            obj.results.strain = permute(strain,[2,3,1]);
            obj.results.stress = permute(stress,[2,3,1]);
            exx = obj.results.strain(:,:,1);
            eyy = obj.results.strain(:,:,2);
            exy = obj.results.strain(:,:,3);
            sxx = obj.results.stress(:,:,1);
            syy = obj.results.stress(:,:,2);
            sxy = obj.results.stress(:,:,3);
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
            obj.results.GPvalues(1,:,:) = exx;
            obj.results.GPvalues(2,:,:) = eyy;
            obj.results.GPvalues(3,:,:) = exy;
            obj.results.GPvalues(4,:,:) = sxx;
            obj.results.GPvalues(5,:,:) = syy;
            obj.results.GPvalues(6,:,:) = sxy;
            obj.results.GPvalues(7,:,:) = e1; 
            obj.results.GPvalues(8,:,:) = e2;
            obj.results.GPvalues(9,:,:) = maxt;
            obj.results.GPvalues(10,:,:) = etheta;
            obj.results.GPvalues(11,:,:) = etr;
            obj.results.GPvalues(12,:,:) = vol;
            obj.results.GPvalues(13,:,:) = s1;
            obj.results.GPvalues(14,:,:) = s2;
            obj.results.GPvalues(15,:,:) = maxs;
            obj.results.GPvalues(16,:,:) = stheta;
            obj.results.GPvalues(17,:,:) = sHM;
            obj.results.GPvalues(18,:,:) = repmat(x,1,nip);
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

