classdef (Abstract) FunctionDual < handle
    methods (Abstract)
        computeValue( x, y );
        computeGradientX( x, y );
        computeGradientY( x, y );
        computeGradientXY( x, y );
    end
    methods
        function [value, gradX, gradY, gradXY] = computeAll( obj, x, y )
            value = computeValue(obj, x, y);
            gradX = computeGradientX(obj,x, y);
            gradY = computeGradientY(obj,x, y);
            gradXY = computeGradientXY(obj,x, y);
        end
    end
end

