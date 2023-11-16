classdef ElastoPlasticity < NonlinearAnalysis 
    
    properties
       
    end
    
    methods
        
        function obj = ElastoPlasticity(felems, mesh)
           obj = obj@NonlinearAnalysis(felems,mesh,"elastoPlasticTangentMatrix");

        end

        function R = computeResidualVector(obj,V)
            R=obj.Pfem(:,1);
            R(:)=0;
            for k=1:max(size(obj.felems))
                %obj.felems{k}.computeStrain( obj.mesh.nodes, obj.qnodal );
                %obj.felems{k}.computeElasticStress()
                %obj.felems{k}.computePlasticStressAndStrainCorrector();
                obj.felems{k}.computeStress( obj.mesh.nodes, obj.qnodal, obj.dq_nodal );
                R = obj.felems{k}.computeInternalForces(mesh.nodes,V,R);
            end
        end
        
        function solve(obj)
            obj.qfem = NewtonRaphsonProcedure(obj);
            obj.qnodal = obj.fromFEMVector(obj.qfem);
        end 

    end
end

