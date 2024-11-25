classdef DOFManagerUniform
    %DOFMANAGERUNIFORM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ndofs;
    end
    
    methods
        function obj = DOFManagerUniform(nnodes,ndofs)
            obj.allDofTypes=ndofs;
            obj.nnodes=nnodes;
        end

        function initDOFs(obj)
        end

        function allocVectors = getAllocationVectors(obj,elems)
            idofs=1:obj.ndofs;
            allocVectors = zeros( size(elems,2)*obj)
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

