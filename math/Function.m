classdef (Abstract) Function < handle
    methods (Abstract)
        computeValue( x );
        computeGradient( x );
    end
    methods
        function [value, grad] = compute( obj, x )
            [ISTILDE, ~, ~, ~] = detectOutputSuppression(nargout);
            if (not(ISTILDE(1))) 
                value = computeValue(obj, x);
            end
            if (not(ISTILDE(2))) 
                grad  = computeGradient(obj,x);
            end
        end
    end
end

