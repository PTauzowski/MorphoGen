classdef PlaneStressElastoPlasticElem < PlaneStressElem
    
    methods
        
        function obj = PlaneStressElastoPlasticElem(sf,elems)
            obj = obj@PlaneStressElem(sf,elems);
            obj.ndofs=["ux" "uy"];
            obj.results.names  = ["exx" "eyy" "exy" "sxx" "syy" "sxy" "e1" "e2" "gmax" "etheta" "tr(e)" "Vol" "s1" "s2" "tmax" "stheta" "sHM" "rho" "pz"];
            obj.results.descriptions  = ["strain member exx" "strain member eyy" "strain member exy" "stress member sxx"...
                 "stress member syy" "stress member sxy" "principal strain e1" "principal strain e2" "maximal shear strain"...
                 "principal strain angle" "strain tensor trace Tr(e)" "volume change" "principal stress s1" "principal stress s2" ...
                 "maximal shear stress" "principal stress angle" "Huber-Mises stress" "Top opt density" "plastic zone"];
        end
  
        function K = elastoPlasticTangentMatrix(obj, nodes, q, varargin)
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
                    Ke = Ke + abs(detJ) * weights(i) * h * B' * obj.mat.tangentD( obj.results.gp.stress(:,k,i), obj.results.gp.dg(k,i) ) * B;
                end
                K(:,:,k) = x(k)*Ke;
            end
            K=K(:);
        end

        function initializeResults(obj)
            nelems = size(obj.elems,1);
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            obj.results.gp.strain = zeros(3,nelems,nip);
            obj.results.gp.pstrain = zeros(3,nelems,nip);
            obj.results.gp.stress = zeros(3,nelems,nip);
            obj.results.gp.dg = zeros(nelems,nip);
            obj.results.gp.all = zeros(size(obj.results.names,2),nelems,nip);
        end
        
        function computeStress(obj, nodes, q, dq, varargin)
            nelems = size(obj.elems,1);
            nnodes = size(obj.elems,2);
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            if ( nargin == 5 ) 
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            B = zeros(3,dim);           
            qelems = reshape( q( obj.elems',:)', nnodes * ndofs, nelems );
            dqelems = reshape( dq( obj.elems',:)', nnodes * ndofs, nelems );
            h = obj.props.h;
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
                    de = B*dqelems(:,k);
                    [ obj.results.gp.stress(:,k,i), obj.results.gp.strain(:,k,i), obj.results.gp.pstrain(:,k,i), obj.results.gp.dg(k,i) ]  = obj.returnMapping( obj.results.gp.strain(:,k,i), obj.results.gp.pstrain(:,k,i), de, x(k)  );
                end
            end            
         end

          function computeResults(obj, nodes, qnodal)
            nelems = size(obj.elems,1);
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);
            exx = obj.results.gp.strain(1,:,:);
            eyy = obj.results.gp.strain(2,:,:);
            exy = obj.results.gp.strain(3,:,:);
            sxx = obj.results.gp.stress(1,:,:);
            syy = obj.results.gp.stress(2,:,:);
            sxy = obj.results.gp.stress(3,:,:);
            e1 = ( exx + eyy ) ./ 2.0 + sqrt( ( (exx - eyy) ./ 2.0 ) .* ( (exx - eyy) ./ 2.0 ) + ( exy .* exy )  );
            e2 = ( exx + eyy ) ./ 2.0 - sqrt( ( (exx - eyy) ./ 2.0 ) .* ( (exx - eyy) ./ 2.0 ) + ( exy .* exy )  );
            maxt = ( e1 - e2 ) ./ 2.0;
            etheta = 2 .* atan( exy ./ (exx - eyy));
            etr = exx + eyy;
            vol = ( 1.0 + e1 ) .* ( 1.0 + e2 );
            s1 = ( sxx + syy ) ./ 2.0 + sqrt( ( (sxx - syy) ./ 2.0 ) .* ( (sxx - syy) ./ 2.0 ) + ( sxy .* sxy )  );
            s2 = ( sxx + syy ) ./ 2.0 - sqrt( ( (sxx - syy) ./ 2.0 ) .* ( (sxx - syy) ./ 2.0 ) + ( sxy .* sxy )  );
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
            obj.results.gp.all(18,:,:) = 1; %repmat(x,1,nip);
            obj.results.gp.all(19,:,:) = obj.results.gp.dg>0;
          end

        function Rfem = computeInternalForces(obj,nodes,refR,V)
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);     
            ndofs = size( obj.ndofs,2);
            nnodes = size(obj.elems,2);
            nelems = size(obj.elems,1);
            dim = nnodes * ndofs;
            V=reshape(V,dim,nelems)';
            Rfem=refR;
            B = zeros(3,dim);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            Rg = zeros( size(B,2), nip );
            h = obj.props.h;
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                %Rfem(V(k,:)) = 0;
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
                    Rfem(V(k,:)) = Rfem(V(k,:)) + abs(detJ) * integrator.weights(i) * h  * B' * obj.results.gp.stress(:,k,i); 
                end
                %Rfem(V(k,:)) = Rfem(V(k,:)) + sum( Rg, 2 );
            end    
        end

        function [ s, e, ep, dg ] = returnMapping( obj, en, epn, de, x  )

            % (i) elastic predictor
            
                E    = x^2*obj.mat.E;
                G    = x^2 * obj.mat.E / 2 / ( 1 + obj.mat.nu );
                eTr  = en  + de; % trial strain
                epTr = epn; % permanent plastic strain
                sTr  = x^2 * obj.mat.D * eTr;
                sy   = x^2 * obj.mat.sy;
            
            % (iI) plasticity condition
                a1 = ( sTr(1) + sTr(2) )^ 2;
                a2 = ( sTr(2) - sTr(1) )^ 2;
                a3 = ( sTr(3) )^ 2;
                
                ksiTr = a1/6 + a2/2 + 2*a3;
                FiTr  = ksiTr/2 - (sy^2) / 3;
                
            %     if  FiTr > 0 && FiTr < 10 
            %         FiTr = 0;
            %     end
                
                if ( FiTr <= 0 )
                    s     = sTr;
                    e     = eTr;
                    ep    = epTr;
                    dg    = 0;
                    ep(:) = 0;
                else
                    [ dg, ksi ] = obj.NR4RetMap( sTr, FiTr, x );
                    A11 = 3 * (1-obj.mat.nu) / (3*(1-obj.mat.nu)+E*dg);
                    A22 = 1 / (1+2*G*dg);
                    A33 = A22;
                    A   = [ 0.5*(A11+A22) 0.5*(A11-A22) 0; 0.5*(A11-A22) 0.5*(A11+A22) 0; 0 0 A33];
                    s   = A * sTr;
                    e   = 1 / x * obj.mat.invD * s;
                    ep  = epTr + dg * sqrt( 2/3 * ksi ); 
                    
                end
        end

        function [ dg, ksi ] = NR4RetMap( obj, sTr, FiTr, x )

            % (1) initial guess
                E  = x^2 * obj.mat.E;
                sy = x^2 * obj.mat.sy;
                nu = obj.mat.nu;
                G  = E / 2 / ( 1 + nu );
                dg  = 0;
                ksi = 0;
                A1 = ( sTr(1) + sTr(2) )^ 2;
                A2 = ( sTr(2) - sTr(1) )^ 2;
                A3 =  sTr(3)^ 2;
            
            % (2) NR iterations
            
                fi = FiTr;
                while( abs(fi/FiTr) >= 0.001 )
                
                      H = 0;
                      ksip = - A1 / (9*(1+E*dg/3/(1-nu))^3) * E / (1-nu)  - 2 * G * ( A2 + 4 * A3 ) / ( 1 + 2 * G * dg )^3;                    
                      Hb = 0;                    
                      Fip = 0.5 * ksip;                    
                      dg = dg - fi / Fip;
                    
                % (3) check for convergence
                      ksi =  A1 / 6 /  (1 + E * dg / 3 / (1-nu) )^2 +  ( 0.5 * A2 + 2 * A3 ) / (1+2*G*dg)^2;
                      fi = 1.0/2.0*ksi - 1.0/3.0 * sy^2;
            
               end
        end
    end
end
