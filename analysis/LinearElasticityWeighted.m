classdef LinearElasticityWeighted < FEAnalysis
   
   properties
        isConst;
   end
   
   methods       
       function obj = LinearElasticityWeighted(felems, mesh, isConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.isConst=isConst;
            obj.rotations=[];
       end
       function K = globalMatrixAggregationWeighted(obj, fname, x)
            K = [];
            ei = getElemIndices(obj);
            for k=1:size(obj.felems,1)
                if ismethod(obj.felems{k},fname)
                    K = [ K obj.felems{k}.(fname)(obj.mesh.nodes,x(ei{k})) ];
                else
                    error("Class " + class(obj.felems{k}) + " or its predecessors not implements function :"+fname);
                end
            end
        end
       function qn = solveWeighted(obj, x)
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           if obj.isConst
                obj.qfem = solver.solve(obj.globalMatrixAggregationWeighted('computeStifnessMatrixConst',x), obj.Pfem);
           else
                obj.qfem = solver.solve(obj.globalMatrixAggregationWeighted('computeStifnessMatrix',x), obj.Pfem);
           end
           obj.qnodal = obj.fromFEMVector(obj.qfem);
           qn=obj.qnodal;
       end
   end
end

