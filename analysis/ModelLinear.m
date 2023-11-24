classdef ModelLinear < handle
   
    properties
        analysis, mesh, fe, x, xp, result_node, result_number;
    end
    
    methods

        function rn = setResultNode(obj, result_node_x)
            rn=obj.mesh.findClosestNode(result_node_x);
            obj.result_node=rn;
        end
        
        function solve(obj)
            obj.analysis.solveWeighted(obj.x);
            obj.analysis.computeElementResults();
        end

        function setupVariables(obj,E,nu,P)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy"], P);
            obj.analysis.solveWeighted(obj.x);
        end

        function u = computeDisplacement(obj,E,nu,P)
            obj.setupVariables(E,nu,P);
            obj.analysis.solveWeighted(obj.x);
            u=obj.analysis.qnodal(obj.result_node,:);
        end

        function sHM = computeHMstress(obj,E,nu,P)
            obj.setupVariables(E,nu,P);
            obj.solve();
            sHM=obj.fe.results.nodal(obj.result_node,obj.result_number);
        end
    end
end

