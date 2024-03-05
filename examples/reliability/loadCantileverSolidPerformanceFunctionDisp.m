classdef  loadCantileverSolidPerformanceFunctionDisp < Function

    properties
        model;
    end

    methods
        function obj = loadCantileverSolidPerformanceFunctionDisp(model)
            obj=obj@Function(3,0.00001)
            obj.model=model;   
        end

        function g = computeValue(obj,points)
            g=0.004+obj.model.computeLinearDisplacement(points);
            % g=zeros(size(points,1),1);
            % for k=1:size(points,1)
            %     ue=obj.model.computeDisplacement(210000,0.3,[points(k,1) points(k,2)]);
            %     g(k)=ue(1);
            % end            
        end

    end

end

