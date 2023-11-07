classdef SinCosApr < Function

     methods
         function obj=SinCosApr()
             obj@Function(2,0.0001);
         end

         function g=computeValue( obj, x )                
                g = sin(x(:,1)).*cos(x(:,2))-0.1;
         end
     end
    
end

