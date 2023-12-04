classdef ElastoPlasticity < NonlinearAnalysis 
    
    properties
       
    end
    
    methods
        
        function obj = ElastoPlasticity(felems, mesh)
           obj = obj@NonlinearAnalysis(felems,mesh,"elastoPlasticTangentMatrix");

        end

        function Fint = computeInternalForces(obj,V)
            refR=obj.Pfem(:,1);
            refR(:)=0;
            cellfun(@(x) x.computeStress( obj.mesh.nodes, obj.qnodal, obj.dq_nodal ), obj.felems);
            F = cellfun( @(x) x.computeInternalForces(obj.mesh.nodes,refR,V), obj.felems,'UniformOutput',false);
            Fint = sum(cell2mat(F),2);
        end
        
        function solve(obj)
            obj.initializeResults();
            NewtonRaphsonProcedure(obj);
            cellfun(@(x) x.computeResults( ),obj.felems);
            obj.qnodal = obj.fromFEMVector(obj.qfem);
            obj.computeElementResults();
        end 

    end
end

