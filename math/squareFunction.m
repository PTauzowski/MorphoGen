classdef squareFunction < Function
    methods
        function y = computeValue( obj, x )
            y=x.^2;
        end
        function dy = computeGradient( obj, x )
            dy=2*x;
        end
    end
end

