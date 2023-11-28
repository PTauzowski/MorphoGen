classdef FEModel
   
    
    properties
        analysis, mesh;
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
            createMaterials();
            applyDirichletBC();
            applyDirichletBC();
        end
       
    end
end

