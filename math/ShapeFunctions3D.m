classdef (Abstract) ShapeFunctions3D < ShapeFunctions
    
    properties
        edges, edgesf;
        faces, fcontours, facesf;
    end

    methods
        function obj = ShapeFunctions()
            obj = obj@ShapeFunctions(3);
        end
    end
    
    methods
        function obj = ShapeFunctions3D()
            obj = obj@ShapeFunctions(3);
        end
    end
    
end

