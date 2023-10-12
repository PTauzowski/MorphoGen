classdef (Abstract) StressIntensityTopologyOptimization < TopologyOptimization
    
    properties
        maxstress, ais, elem_list, maxais, penal, xr, xfat, pnormfat;
    end

    methods
            
        function obj = StressIntensityTopologyOptimization(numberOfConstraints,Rmin,FEproblem,maxais,penal,is_const)
            obj=obj@TopologyOptimization(numberOfConstraints,Rmin,FEproblem,is_const);
            obj.maxais=maxais;
            obj.penal=penal;
            obj.x(:)=1;
            obj.xfat=obj.x;
            obj.pnormfat=100;
        end

        function ret = isNotFinished(obj)
            obj.constrValues = obj.computeInequalityConstraints(obj.x);
            if find( obj.constrValues < 0 , 1 )
                    ret = false;
                    return;
            end
            ret=true;
        end

        function updateDesign(obj)
            obj.ais = obj.weights * obj.computeAverageIntensities();
            obj.removeStressed();
            wais=obj.ais/max(obj.ais);
            obj.x(obj.erased_elems) = min(1, max( 0.001, obj.x(obj.erased_elems) .* wais(obj.erased_elems) .^ obj.penal)); 
            obj.x( not(obj.erased_elems) ) = 1; 
            obj.x( obj.const_elems ) = 1;
        end

        function removeStressed( obj )
            notErasedID = find( not( obj.erased_elems )  );
            notErasedID = setxor(notErasedID,intersect(notErasedID, obj.const_elems));

            maxaisprc =(max(obj.ais(notErasedID)) - min(obj.ais(notErasedID)) ) * obj.maxais;
            obj.elem_list = obj.ais(notErasedID) < min(obj.ais(notErasedID)) + maxaisprc;

            obj.elem_list = notErasedID( obj.elem_list );
            obj.erased_elems( obj.elem_list ) = true;
        end

        function ais = computeAverageIntensities(obj)
            obj.qnodal = obj.FEproblem.solveWeighted((obj.x).^obj.penal);
            obj.FEproblem.computeElementResults(obj.qnodal,obj.x.^obj.penal);
            ais = zeros(obj.FEproblem.getTotalElemsNumber(),1);
            for i=1:size(obj.FEproblem.felems,1)
               hmIndex=find(obj.FEproblem.felems{i}.results.names == "sHM");
               for j=1:size(obj.FEproblem.felems{i}.elems,1)
                    %ais(obj.elem_inds{i}(j)) = mean( obj.linearElasticProblem.felems{i}.results.GPvalues(hmIndex,j,:) );
                    ais(obj.elem_inds{i}(j)) = mean( obj.FEproblem.felems{i}.results.nodal(obj.FEproblem.felems{i}.elems(j,:),hmIndex) );
               end
            end
            obj.maxstress = [ obj.maxstress max(ais) ];
            ais = ais / max(ais);
        end

    end
end

