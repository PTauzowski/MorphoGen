classdef (Abstract) StressIntensityTopologyOptimization < TopologyOptimization
    
    properties
        maxstress, ais, elem_list, maxais, penal, xr, xfat, pnormfat, allx;
    end

    methods
            
        function obj = StressIntensityTopologyOptimization(numberOfConstraints,Rmin,FEAnalysis,maxais,penal,is_const)
            obj=obj@TopologyOptimization(numberOfConstraints,Rmin,FEAnalysis,is_const);
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

        function resetAnalysis(obj)
            obj.x(:)=1;
            obj.totalFENumber = obj.FEAnalysis.getTotalElemsNumber();
            obj.erased_elems = false(obj.totalFENumber,1);
        end

        function updateDesign(obj)
            obj.ais = obj.weights * obj.computeAverageIntensities();
            obj.removeStressed();
            wais=obj.ais/max(obj.ais);
            obj.x(obj.erased_elems) = min(1, max( 0.001, obj.x(obj.erased_elems) .* wais(obj.erased_elems) .^ obj.penal)); 
            obj.x( not(obj.erased_elems) ) = 1; 
            obj.x( obj.const_elems ) = 1;
            obj.allx=[obj.allx obj.x];
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
            obj.qnodal = obj.FEAnalysis.solveWeighted((obj.x).^obj.penal);
            obj.FEAnalysis.computeElementResults(obj.x.^obj.penal);
            ais = zeros(obj.FEAnalysis.getTotalElemsNumber(),1);
            for i=1:size(obj.FEAnalysis.felems,1)
               hmIndex=find(obj.FEAnalysis.felems{i}.results.names == "sHM");
               for j=1:size(obj.FEAnalysis.felems{i}.elems,1)
                    %ais(obj.elem_inds{i}(j)) = mean( obj.linearElasticProblem.felems{i}.results.GPvalues(hmIndex,j,:) );
                    ais(obj.elem_inds{i}(j)) = mean( obj.FEAnalysis.felems{i}.results.nodal.all(obj.FEAnalysis.felems{i}.elems(j,:),hmIndex) );
               end
            end
            obj.maxstress = [ obj.maxstress max(ais) ];
            ais = ais / max(ais);
        end

        function setFrame( obj, k )
            if k > 0 && k <= size( obj.allx,2 )
                obj.x=obj.allx(:,k);
                obj.FEAnalysis.computeElementResults(obj.x.^obj.penal);
                return;
            end
            fprintf("The specified iteration number %d is outside the allowed range %d - %d\n",k,1,size( obj.allx,2 ));
        end

        function ifound = findFrame(obj, vol )
            i1=1; 
            i2=size(obj.allx,2);
            v1 = sum(obj.allx(:,1))/size(obj.allx,1);
            v2 = sum(obj.allx(:,end))/size(obj.allx,1);
            if vol > v1 || vol < v2
                fprintf("The specified volume fraction %1.3f is outside the existing range %1.3f - %1.3f\n",vol,v1,v2);
            else
                while true
                    i=round((i1+i2)/2);
                    v=sum(obj.allx(:,i))/size(obj.allx,1);
                    %fprintf("Mid point i=%d, v=%1.3f\n",i,v);
                    if i2-i1==1
                        ifound=i1;
                        break;
                    end
                    if v>vol
                        v1=v;
                        i1=i;
                    else
                        v2=v;
                        i2=i;
                    end
                end
            end
            
        end

        function printVolumeFractionRanges(obj)
            i1=1; 
            i2=size(obj.allx,2);
            obj.setFrame(i1);
            [v1, va1, vc1] = obj.computeVolumeFraction();
            obj.setFrame(i2);
            [v2, va2, vc2] = obj.computeVolumeFraction();
            fprintf("Frame %5d, VolFr=%1.3f, VolFr_active=%1.3f, VolFr_const=%1.3f\n",i1,v1,va1,vc1);
            fprintf("Frame %5d, VolFr=%1.3f, VolFr_active=%1.3f, VolFr_const=%1.3f\n",i2,v2,va2,vc2);
        end

    end
end

