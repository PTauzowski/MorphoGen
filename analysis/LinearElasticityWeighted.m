classdef LinearElasticityWeighted < FEAnalysis
   
   properties
        isConst,selfLoadFactor;
   end
   
   methods       
       function obj = LinearElasticityWeighted(felems, mesh, isConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.isConst=isConst;
            obj.rotations=[];
            obj.selfLoadFactor=-1;
       end
       % function K = globalMatrixAggregationWeighted(obj, fname, x)
       %      K = [];
       %      ei = getElemIndices(obj);
       %      for k=1:size(obj.felems,2)
       %          if ismethod(obj.felems{k},fname)
       %              K = [ K; obj.felems{k}.(fname)(obj.mesh.nodes,x(ei{k})), [] ];
       %          else
       %              error("Class " + class(obj.felems{k}) + " or its predecessors not implements function :"+fname);
       %          end
       %      end
       %  end
       function qfem = solve(obj, x)
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
           if obj.selfLoadFactor>0
               obj.Pnodal(:) = 0;
               obj.loadElementsSelfWeight(x,0.1);
               obj.Pfem = obj.selfLoadFactor .* obj.toFEMVector(obj.Pnodal);
               obj.Pnodal(:) = 0;
           end
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
            end
            
            obj.qfem = solver.solve(obj.assemblyGlobalMatrix('computeStifnessMatrix',x,obj.isConst), obj.Pfem);

           qfem=obj.qfem;
           obj.qnodal=obj.fromFEMVector(qfem(:,1));
       end
   end
end

