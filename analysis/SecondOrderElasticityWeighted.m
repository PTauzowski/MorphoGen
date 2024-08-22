classdef SecondOrderElasticityWeighted < FEAnalysis
   
   properties
        isConst,lambda;
   end
   
   methods       
       function obj = SecondOrderElasticityWeighted(felems, mesh, lambda, isConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.lambda = lambda;
            obj.isConst=isConst;
            obj.rotations=[];
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
           obj.qfem = solver.solve(K+obj.lambda*Kg, obj.Pfem);
           obj.qnodal=obj.fromFEMVector(obj.qfem(:,1));
           qfem=obj.qfem;
       end

       
   end
end

