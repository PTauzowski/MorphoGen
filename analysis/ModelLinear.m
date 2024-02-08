classdef ModelLinear < FEModel
   
    properties
        analysis, fe, xp, result_node, result_number;
    end

    methods

        function setX(obj,x)
            obj.x=x;
        end

         function x=getX(obj)
            x=obj.x;
        end

        function rn = setResultNode(obj, result_node_x)
            rn=obj.mesh.findClosestNode(result_node_x);
            obj.result_node=rn;
        end
        
        function solveWeighted(obj)
            obj.analysis.solveWeighted(obj.x);
            obj.analysis.initializeResults();
            obj.analysis.computeElementResults();
        end
    
        function setupLoad(obj,P)
            obj.analysis.clearCurrentLoad();
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], P);
        end
    
        function setupVariables(obj,E,nu,P)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], P);
        end
    
        function u = computeDisplacement(obj,E,nu,pressure)
            obj.analysis.clearCurrentLoad();
            obj.setupVariables(E,nu,pressure);
            obj.analysis.solveWeighted(obj.x);
            u=obj.analysis.qnodal(obj.result_node,:);
        end
    
        function sHM = computeHMstress(obj,E,nu,pressure)
            obj.analysis.clearCurrentLoad();
            obj.setupVariables(E,nu,pressure);
            obj.analysis.solveWeighted(obj.x);
            obj.analysis.computeElementResults(obj.x);
            sHM=obj.fe.results.nodal.all(obj.result_node,obj.result_number);
        end 
    
        function pstress = computePenalizedHMstress(obj,E,nu,pressure,penalty)
            obj.analysis.clearCurrentLoad();
            obj.setupVariables(E,nu,pressure);
            obj.analysis.solveWeighted(obj.x);
            obj.analysis.computeElementResults(obj.x);
            pstress=sum(obj.fe.results.nodal.all(:,18).*obj.fe.results.nodal.all(:,obj.result_number))^(1/penalty);
            %pstress=sum(obj.fe.results.nodal.all(:,obj.result_number))^(1/penalty);
        end
    
        function plotModel( obj )
            obj.analysis.plotFiniteElements();
            %obj.analysis.plotCurrentLoad();
            %obj.analysis.plotSupport();
            plot(obj.mesh.nodes(obj.result_node,1),obj.mesh.nodes(obj.result_node,2),"Marker","o");
        end
    
        function plotNodes(obj)
            plot(obj.mesh.nodes(:,1),obj.mesh.nodes(:,2),'.');
        end
    
        function vol = getVolumeFraction(obj)
            vol=sum(obj.x(:))/size(obj.x,1);
        end

    end
end

