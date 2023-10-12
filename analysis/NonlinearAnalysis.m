classdef NonlinearAnalysis < FEAnalysis
    
    properties
        convEnd, iteration;
        tangentMatrixProcedureName;
    end
    
    methods (Abstract)
        computeResidualVector();
    end
    
    methods
        function obj = NonlinearAnalysis(felems,mesh,tangentMatrixProcedureName)
           obj = obj@FEAnalysis(felems,mesh);
           obj.tangentMatrixProcedureName=tangentMatrixProcedureName;
           obj.convEnd = 1.0E-06;
        end
        
        function qfem = NewtonRaphsonProcedure(obj)
            conv=1.0;
            [I,J,~,~] = obj.globalMatrixIndices();
            solver =  LinearEquationsSystem(I,J,obj.supports);
            obj.prepareRHSVectors();
            obj.iteraton=1;
            dP = obj.Pfem(:,0);
            while conv > obj.convEnd
                Kt = obj.globalSolutionDependendMatrixAggregation( obj.tangentMatrixProcedureName );
                obj.qfem = solver.solve( Kt, dP );
                obj.qnodal = obj.fromFEMVector( obj.qfem );
                obj.computeElementResults( obj.qnodal );
                R=obj.computeResidualVector();
                dP=obj.Pfem(:,0)-R;
                conv  = norm( dP ) / norm( R );
                obj.iteration = obj.iteration + 1;
            end
        end
    end
end

