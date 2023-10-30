classdef  matPerformanceFunction < Function

    properties
        analysis, material, nelem, nnode, nres;
    end

    methods
        function obj = matPerformanceFunction(analysis,material,nnode,nelem,nres)
            obj=obj@Function(2,[20 0.0001]);
            obj.analysis = analysis;
            obj.material=material;
            obj.nnode=nnode;
            obj.nelem=nelem;
            obj.nres=nres;
        end

        function g = computeValue(obj,x)
            g=zeros(size(x,1));
            for k=1:size(x,1)
                obj.material.setElasticIzo(x(1), x(2));
                obj.analysis.solve();
               % g(k)=obj.analysis.felems{obj.nelem}.results.nodal(obj.nnode,obj.nres);
                g(k)=obj.analysis.qnodal(obj.nnode,2);
            end
        end

    end

end

