classdef ElasticHarmonicVibrations < FEAnalysis
   
   properties
        isMeshConst, isLoadConst, lambdas, lambdas1, omegas, frequencies, mode, modes, modes1, nmodes, count, P0, Mnodal, correlation_matrix, freedofs, loadVectors,save_frame_modes,M1;
   end
   
   methods       
       
       function obj = ElasticHarmonicVibrations(felems, mesh, mode, isMeshConst, isLoadConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.mode = mode;
            obj.isMeshConst=isMeshConst;
            obj.isLoadConst=isLoadConst;
            obj.rotations=[];
            obj.count=0;
            obj.modes=[];
            obj.lambdas=[];
            obj.omegas=[];
            obj.frequencies=[];
            obj.loadVectors=[];
            obj.nmodes=30;
            obj.save_frame_modes=false;
            obj.correlation_matrix=zeros(obj.nmodes,obj.nmodes);
            obj.Mnodal = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
       end

       function massClosestNode(obj, x, dofnames, values )
           m_nodes = obj.mesh.findClosestNode(x);
           obj.Mnodal( m_nodes, obj.findDOFsIndices( dofnames ) ) = obj.Mnodal( m_nodes, obj.findDOFsIndices( dofnames ) ) + values;
        end

       function [K, M] = computeMatrices(obj,x)
           K  = obj.createSparseMatrix( obj.assemblyGlobalMatrix('computeStifnessMatrix', x, obj.isMeshConst) );
           M  = obj.createSparseMatrix( obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isMeshConst) );
       end

       function vM = addMassLumped(obj, vM, I, J)
           Mfem = obj.toFEMVector(obj.Mnodal);
           diag_idx = I==J;
           lumped_mass_dof = find(Mfem );
           for i=1:numel(lumped_mass_dof)
               loaded_idx = find(I(diag_idx)==lumped_mass_dof(i));
               vM(loaded_idx) = vM(loaded_idx) + vM(loaded_idx).*Mfem(I(loaded_idx));
           end
       end

       function qfem = solve(obj, x)
           obj.count=obj.count+1;
           [I,J,~] = obj.globalMatrixIndices();
%           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           vK  = obj.assemblyGlobalMatrix('computeStifnessMatrix', x, obj.isMeshConst);
           vM  = obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isMeshConst);
           vM = obj.addMassLumped(vM,I,J);
           obj.freedofs=solver.freedofs;

            K = solver.createSparseMatrix(vK);
            M = solver.createSparseMatrix(vM);
           % save("LoadForce.mat",  "K", "M", "x");

           [modes, lambdas] = solver.solveEigenproblem(vK,vM,obj.nmodes);
           lambdaVec = diag(lambdas);
            omegaVec = sqrt(lambdaVec);
            freqVec = omegaVec / (2*pi);
            
            if obj.save_frame_modes
                % Store history of modes: big memory cost
                obj.modes       = cat(3, obj.modes, modes);
                obj.lambdas     = [obj.lambdas lambdaVec];
                obj.omegas      = [obj.omegas omegaVec];
                obj.frequencies = [obj.frequencies freqVec];
            else
                % Keep only current modes: memory-friendly
                obj.modes       = modes;
                obj.lambdas     = lambdaVec;
                obj.omegas      = omegaVec;
                obj.frequencies = freqVec;
            end
           if obj.count==1
               obj.modes1=modes;
               obj.lambdas1=lambdas;
               obj.M1=M;
           end

           obj.Pfem=zeros(solver.dim,1);
           obj.qfem=zeros(solver.dim,1);

           if obj.isLoadConst
                obj.Pfem(solver.freedofs) =  obj.lambdas1(obj.mode,obj.mode) * obj.M1 * obj.modes1(solver.freedofs,obj.mode);
           else
                obj.Pfem(solver.freedofs) =  lambdas(obj.mode,obj.mode)*M*modes(solver.freedofs,obj.mode);         
           end
           obj.qfem = solver.solve(vK, obj.Pfem);
           obj.loadVectors = [ obj.loadVectors obj.Pfem ];
           obj.qnodal=obj.fromFEMVector(obj.qfem(:,1));
           qfem=obj.qfem;

       end

       function Pfem = computeLoadVector(obj, frame, x ) 
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           vM  = obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isMeshConst);
           M = solver.createSparseMatrix(vM);
           Pfem=0*obj.Pfem;
           Pfem(solver.freedofs) = obj.lambdas(obj.mode,1)*M*obj.modes(solver.freedofs,obj.mode,1);
       end

       function Pfem = getLoadVector(obj, frame ) 
           Pfem = obj.loadVectors(:, frame);
       end

       function conds = tabMatrixCondition( obj, alphas, mode )
           conds=[];
           disp(['Matrix conditions for mode: ' num2str(mode)]);
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           K  =  obj.assemblyGlobalMatrix('computeStifnessMatrix', 1, obj.isConst);
           M  =  obj.assemblyGlobalMatrix('computeMassMatrix', 1, obj.isConst);
           [~, lambdas] = solver.solveEigenproblem(K,M,10);
           lambdas = diag(lambdas(1:5,1:5));

           for k=1:numel(alphas)
               disp(['iter: ' num2str(k)]);
               conds = [conds rcond(full(solver.createSparseMatrix(K-alphas(k)*lambdas(mode)*M)))];
           end
       end
        
       function computeSelfCorrelationMatrix(obj)
           for k=1:obj.nmodes
               for l=1:obj.nmodes
                    obj.correlation_matrix(k,l)=abs(obj.modes(:,k)'*obj.modes(:,l))/norm(obj.modes(:,k))/norm(obj.modes(:,l));
               end
           end
       end

       function Corr = computeCorrelationMatrix(obj, dd_modes, nmodes)
           Corr=zeros(nmodes,nmodes);
           for k=1:nmodes
               for l=1:nmodes
                    Corr(k,l)=abs(dd_modes(:,k)'*obj.modes(:,l,end))/norm(dd_modes(:,k))/norm(obj.modes(:,l,end));
               end
           end
       end

       function cv = computeCorrelationVector(obj,test_mode)
           cv=zeros(1,obj.nmodes);
           for k=1:obj.nmodes
                    cv(k)=abs(test_mode'*obj.modes(:,k,end))/norm(test_mode)/norm(obj.modes(:,k,end));
           end
       end

   end
end

