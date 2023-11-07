classdef SpatialObject < Function
   
    properties
        ndiv
    end

    methods (Abstract)
        paint(obj);
    end 
    
    methods
        
        function obj = SpatialObject(dim,ndiv)
            obj=obj@Function(dim,0.001);
            obj.ndiv = ndiv;
        end    

    end
end

