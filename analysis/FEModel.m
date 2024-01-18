classdef FEModel
  
    properties
       mesh, x, rho, supports, rotations;
    end

    methods (Abstract)
        generateMesh();
        createFiniteElements();
        createMaterials();
        applyDirichletBC();
        applyNeumannBC();
    end
    
    methods
        function obj = FEModel()
            obj.generateMesh();
            obj.createFiniteElements();
            obj.createMaterials();
            obj.applyDirichletBC();
            obj.applyDirichletBC();
        end
       
    end
end

