classdef Selector
        
    properties
        fn;
        tolerance;
    end
    
    methods
        function obj = Selector(fn,varargin)
            obj.fn=fn;
            obj.tolerance=1.0E-04;
            if nargin==2
                obj.tolerance = varargin{1};
            end
        end
        
        function s = select( obj, points )
            if islogical(obj.fn)
                s=obj.fn;
            else
                s = abs(obj.fn(points)) < obj.tolerance;
            end
           
        end
        
    end
end

