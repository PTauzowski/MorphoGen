classdef SecondOrderDynamicElasticityWeighted < FEAnalysis
   
   properties
        isConst,lambda, omegas, lambdac, modes1, modes2, modes3, modes4, modes5;
   end
   
   methods       
       function obj = SecondOrderDynamicElasticityWeighted(felems, mesh, lambdac, isConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.lambdac = lambdac;
            obj.isConst=isConst;
            obj.rotations=[];
            obj.modes1=[];
            obj.modes2=[];
            obj.modes3=[];
            obj.modes4=[];
            obj.modes5=[];

       end
       function K = globalMatrixAggregationWeighted(obj, fname, x)
            K = [];
            ei = getElemIndices(obj);
            for k=1:size(obj.felems,2)
                if ismethod(obj.felems{k},fname)
                    K = [ K; obj.felems{k}.(fname)(obj.mesh.nodes,x(ei{k})) ];
                else
                    error("Class " + class(obj.felems{k}) + " or its predecessors not implements function :"+fname);
                end
            end
        end
       function qfem = solveWeighted(obj, x)
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           if obj.isConst
               K = obj.globalMatrixAggregationWeighted('computeStifnessMatrixConst',x);
           else
               K = obj.globalMatrixAggregationWeighted('computeStifnessMatrix',x);
           end
           obj.qfem = solver.solve(K, obj.Pfem);
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults(x);
           Kg = obj.globalMatrixAggregationWeighted('computeGeometricStifnessMatrix',x);
           M  = obj.globalMatrixAggregationWeighted('computeMassMatrix',x);
           [eigenvectors, lambdas] = solver.solveEigenproblem(K,Kg,10);
           [mode, omegas] = solver.solveEigenproblem(K,M,10);
           obj.modes1=[ obj.modes1 mode(:,1) ];
           obj.modes2=[ obj.modes2 mode(:,2) ];
           obj.modes3=[ obj.modes3 mode(:,3) ];
           obj.modes4=[ obj.modes4 mode(:,4) ];
           obj.modes5=[ obj.modes5 mode(:,5) ];
           obj.lambda=lambdas(1);
           obj.omegas=diag(omegas(1:5,1:5));
           obj.qfem=0*obj.Pfem;
           %obj.qfem(solver.freedofs)=eigenvectors(:,1);
           obj.qfem = solver.solve(K-obj.lambdac*abs(lambdas(1))*Kg, obj.Pfem);
           obj.qnodal=obj.fromFEMVector(obj.qfem(:,1));
           qfem=obj.qfem;
       end

       function [plCor] = getModesCorrelation(obj,mode)
           plCor=[];
           for k=2:size(mode,2)
               plCor = [plCor abs(mode(:,k-1)'*mode(:,k))/norm(mode(:,k-1))/norm(mode(:,k))];
           end
       end

       
   end
end

