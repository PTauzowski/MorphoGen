classdef (Abstract) ShapeFunctions < Function
    
    properties
          localNodes;
          pattern;
          vertices;
    end
    
    methods(Abstract)
        createIntegrator(varargin);
    end
    
    methods
        function obj = ShapeFunctions(dim)
            obj = obj@Function(dim,0.00001);
        end
        function isCorrect = selfTest( obj )
            isCorrect = isdiag( computeValue( obj.localNodes ) );
        end
        function N = getRecoveryMatrix(obj)
            integrator = obj.createIntegrator();
            N = obj.computeValue( integrator.points );
        end
    end
end

