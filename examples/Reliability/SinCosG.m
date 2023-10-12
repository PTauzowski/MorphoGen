classdef SinCosG < Function

     methods
         function g=computeValue( obj, x )
                g = sin(x(:,1)).*cos(x(:,2));
         end
         function grad = computeGradient( obj, x )
             grad = [cos(x(:,1)).*cos(x(:,2)) -sin(x(:,1)).*sin(x(:,2))];
         end
     end
    
end

