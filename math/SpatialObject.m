classdef SpatialObject < Function
   
    properties
        ndiv
    end

    methods (Abstract)
        paint(obj);
    end 
    
    methods
        
        function obj = SpatialObject(ndiv)
            obj.ndiv = ndiv;
        end    

    end
end

