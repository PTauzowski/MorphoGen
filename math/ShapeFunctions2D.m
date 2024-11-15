classdef (Abstract) ShapeFunctions2D < ShapeFunctions
    
    properties
        edges, edgesf, contour;
    end

    methods
        function obj = ShapeFunctions2D()
            obj = obj@ShapeFunctions(2);
        end
    end

end

