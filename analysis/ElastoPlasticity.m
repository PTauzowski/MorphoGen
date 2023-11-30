classdef ElastoPlasticity < NonlinearAnalysis 
    
    properties
       
    end
    
    methods
        
        function obj = ElastoPlasticity(felems, mesh)
           obj = obj@NonlinearAnalysis(felems,mesh,"elastoPlasticTangentMatrix");

        end

        function R = computeResidualVector(obj,V)
            refR=obj.Pfem(:,1);
            refR(:)=0;
            % for k=1:max(size(obj.felems))
            %     obj.felems{k}.computeResults( obj.mesh.nodes, obj.qnodal, obj.dq_nodal );
            %     R = obj.felems{k}.computeInternalForces(mesh.nodes,V,R);
            % end
            cellfun(@(x) x.computeStress( obj.mesh.nodes, obj.qnodal, obj.dq_nodal ), obj.felems);
            R = sum(cellfun( @(x) x.computeInternalForces(obj.mesh.nodes,V,refR), obj.felems));
        end
        
        function solve(obj)
            obj.initializeResults();
            obj.qfem = NewtonRaphsonProcedure(obj);
            cellfun(@(x) x.computeResults( ),obj.felems);
            obj.qnodal = obj.fromFEMVector(obj.qfem);
        end 

    end
end

