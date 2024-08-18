classdef Frame3D < FiniteElement
      
    methods
        function obj = Frame3D(elems)
            obj = obj@FiniteElement(ShapeFunctionsFrame3D,elems);
        end
    end
end


