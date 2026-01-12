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
            if numel(x) == 1
                x = ones(obj.getTotalElemsNumber() ,1);
            end
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
                obj.qfem = solver.solve(obj.globalMatrixAggregationWeighted('computeStifnessMatrixConst',x), obj.Pfem);
           else
                obj.qfem = solver.solve(obj.globalMatrixAggregationWeighted('computeStifnessMatrix',x), obj.Pfem);
           end
           qfem=obj.qfem;
           obj.qnodal=obj.fromFEMVector(qfem(:,1));
       end

       function obj = saveMatrices(obj,filename)
          [I,J,~] = obj.globalMatrixIndices();
          if obj.isConst
                K = sparse(I,J,obj.globalMatrixAggregationWeighted('computeStifnessMatrixConst',1));
          else
                K = sparse(I,J,obj.globalMatrixAggregationWeighted('computeStifnessMatrix',1));
          end
           supports = obj.supports;
           obj.prepareRHSVectors();
           P = obj.Pfem;
           nodes = obj.mesh.nodes;
           elems = obj.mesh.elems;
           save(filename, "nodes","elems","K","P","supports");
       end
   end
end

