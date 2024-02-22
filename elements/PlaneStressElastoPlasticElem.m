classdef PlaneStressElastoPlasticElem < PlaneStressElem
    
    
    properties
        Property1
    end
    
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

         function computeResults(obj)
            nelems = size(obj.elems,1);
            nip = size(integrator.points,1);
            obj.results.GPvalues = zeros(size(obj.results.names,2),nelems,nip);
            exx = obj.results.strain(:,:,1);
            eyy = obj.results.strain(:,:,2);
            exy = obj.results.strain(:,:,3);
            sxx = obj.results.stress(:,:,1);
            syy = obj.results.stress(:,:,2);
            sxy = obj.results.stress(:,:,3);
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
            obj.results.GPvalues(19,:,:) = obj.resultsdg>0;
         end

         function computeStrain(obj, nodes, q)
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
            B = zeros(3,dim);
            strain = zeros(3,nelems,nip);
            strainp = zeros(3,nelems,nip);
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
                    e = reshape( q(obj.elems(k,:),:)', nnodes * ndofs,1 );
                    s = x(k)*D*e;
                    strain(:,k,i) = B*qelems(:,k);
                end
            end
            obj.results.strain = permute(strain,[2,3,1]);
            obj.results.strainp = permute(strain,[2,3,1]);
            obj.results.strainp(:) = 0;
            obj.results.dg=zeros(nelems,nip);
            
         end

        function computeElasticStress(obj, varargin)
            nelems = size(obj.elems,1);
            if ( nargin == 4 )
                x=varargin{1};
            else
                x=ones(nelems,1);
            end
            obj.results.stress = zeros(nelems,nip,3);
            D = obj.mat.D;
            for k=1:nelems
                for i=1:nip
                    obj.results.stress(:,k,i) = x(k)*D*obj.results.strain(k,i,:);
                end
            end
        end

        function computePlasticStressAndStrainCorrector(obj)
            integrator = obj.sf.createIntegrator();
            nelems = size(obj.elems,1);
            nip = size(integrator.points,1);     
            for k=1:nelems
                for i=1:nip
                    [ obj.results.GPvalues.stress(k,i,:), obj.results.GPvalues.strain(k,i,:), obj.results.GPvalues.strainp(k,i,:), obj.results.GPvalues.dg(k,i) ] = elemClass.returnMapping( obj.results.GPvalues.strain(k,i,:), obj.results.GPvalues.strainp(k,i,:), de(:,i), mat, x );
                end
            end
        end

        function Rfem = computeInternalForces(obj,nodes,Rfem,V)
            integrator = obj.sf.createIntegrator();
            nip = size(integrator.points,1);     
            ndofs = size( obj.ndofs,2);
            dim = nnodes * ndofs;
            B = zeros(3,dim);
            dN = obj.sf.computeGradient( integrator.points );
            nnd = size(dN,1); 
            dNtr = permute(dN,[2,1,3]);
            dNtrc = cell(size(dNtr,3),1);
            for i=1:nip
                dNtrc{i}=dNtr(:,:,i);
            end
            Rg = zeros( size(B,2), ngp );
            h = obj.props.h;
            for k=1:nelems
                elemX = nodes(obj.elems(k,:),:);
                Rfem(V(k,:)) = 0;
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
                    Rg(:,i)= abs(detJ) * integrator.weights(i) * h  * B' * obj.results.GPvalues.stress(k,i,:); 
                end
                Rfem(V(k,:)) = Rfem(V(k,:)) + sum( Rg, 2 );
            end    
        end

        function [ s, e, ep, dg ] = returnMapping( en, epn, de, x  )

            % (i) elastic predictor
            
                G    = x * mat.E / 2 / ( 1 + mat.nu );
                eTr  = en  + de; % trial strain
                epTr = epn; % permanent plastic strain
                sTr  = x * mat.D * eTr;
                sy   = mat.sy;
            
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
                    [ dg, ksi ] = NR4RetMap( sTr, FiTr, mat, x );
                    A11 = 3 * (1-mat.nu) / (3*(1-mat.nu)+x*mat.E*dg);
                    A22 = 1 / (1+2*G*dg);
                    A33 = A22;
                    A   = [ 0.5*(A11+A22) 0.5*(A11-A22) 0; 0.5*(A11-A22) 0.5*(A11+A22) 0; 0 0 A33];
                    s   = A * sTr;
                    e   = 1 / x * mat.invD * s;
                    ep  = epTr + dg * sqrt( 2/3 * ksi ); 
                    
                end
        end

        function [ dg, ksi ] = NR4RetMap( sTr, FiTr, mat, x )

            % (1) initial guess
                E  = x * mat.E;
                sy = mat.sy;
                nu = mat.nu;
                G  = E / 2 / ( 1 + nu );
                dg  = 0;
                ksi = 0;
                A1 = ( sTr(1) + sTr(2) )^ 2;
                A2 = ( sTr(2) - sTr(1) )^ 2;
                A3 =  sTr(3)^ 2;
            
            % (2) NR iterations
            
                fi = FiTr;
                while( abs(fi/FiTr) >= 0.5 )
                
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

