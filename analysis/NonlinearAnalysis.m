classdef NonlinearAnalysis < FEAnalysis
    
    properties
        convEnd, iteration, dq_fem, dq_nodal;
        tangentMatrixProcedureName;
    end
    
    methods (Abstract)
        computeInternalForces();
    end
    
    methods
        function obj = NonlinearAnalysis(felems,mesh,tangentMatrixProcedureName)
           obj = obj@FEAnalysis(felems,mesh);
           obj.tangentMatrixProcedureName=tangentMatrixProcedureName;
           obj.convEnd = 1.0E-06;
        end
        
        function qfem = NewtonRaphsonProcedure(obj)
            conv=1.0;
            [I,J,V,~] = obj.globalMatrixIndices();
            solver =  LinearEquationsSystem(I,J,obj.supports);
            obj.prepareRHSVectors();
            obj.iteration=1;
            dP = obj.getCurrentFEMlLoad();
            P=obj.Pfem(:,1);
            while conv > obj.convEnd
                Kt = obj.globalSolutionDependendMatrixAggregation( obj.tangentMatrixProcedureName );
                obj.dq_fem = solver.solveClassical( Kt, dP );
                if obj.iteration==1
                    obj.qfem = obj.dq_fem;
                else
                    obj.qfem = obj.qfem + obj.dq_fem;
                end
                obj.qnodal = obj.fromFEMVector( obj.qfem );
                obj.dq_nodal = obj.fromFEMVector( obj.dq_fem );
                Fint = obj.computeInternalForces(V);
                R = P-Fint;
                conv = norm( R(solver.freedofs ) ) / norm( P( solver.freedofs ) )
                obj.iteration = obj.iteration + 1;
            end
        end
    end
end

