classdef  loadPerformanceFunction < Function

    properties
        analysis, material, nelem, nnode, nres,loadedEdgeSelectorX,loadedEdgeSelectorY;
    end

    methods
        function obj = loadPerformanceFunction(analysis,material,nnode,nelem,nres)
            obj=obj@Function(2,[0.001 0.001]);
            obj.analysis = analysis;
            obj.material=material;
            obj.nnode=nnode;
            obj.nelem=nelem;
            obj.nres=nres;
        end

        function g = computeValue(obj,points)
            g=zeros(size(points,1),1);
            for k=1:size(points,1)
                obj.analysis.clearCurrentLoad();
                obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdgeSelectorX, ["ux" "uy"], @(x)( x*0 + [ points(k,1) 0 ] ));
                obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdgeSelectorY, ["ux" "uy"], @(x)( x*0 + [ 0 points(k,2) ] ));
                obj.analysis.solve();
                %g(k)=obj.analysis.felems{obj.nelem}.results.nodal(obj.nnode,obj.nres);
                g(k)=obj.analysis.qnodal(obj.nnode,2);
            end            
        end

    end

end

