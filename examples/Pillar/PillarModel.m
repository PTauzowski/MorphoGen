classdef PillarModel < ModelLinear
    

    properties
        top_R
        pillar_layers, ground_layers
        pillar_res, ground_res
        pillar_chem, ground_chem
        int_th
        sf
    end

    methods
        function obj = PillarModel( top_R, pillar_layers, ground_layers, pillar_res, ground_res, pillar_chem, ground_chem, int_th, sf )
            obj.top_R = top_R;
            obj.pillar_layers = pillar_layers;
            obj.ground_layers = ground_layers;
            obj.pillar_res = pillar_res;
            obj.ground_res = ground_res;
            obj.pillar_chem = pillar_chem;
            obj.ground_chem = ground_chem;
            obj.int_th = int_th;
            obj.sf = sf;
            obj.generateMesh();
           
        end

        function obj = generateMesh(obj)
            obj.mesh = Mesh();
            pillar_height = sum(obj.pillar_layers);
            ground_depth  = sum(obj.ground_layers);
            nrXY = 6;
        
            Rbase = obj.top_R + pillar_height * tan(5*pi/180);
            Rtop  = obj.top_R;
        
            % ground matches pillar base
            obj.mesh.addLayeredQuarterCylinder([0, 0, -ground_depth], Rbase, obj.ground_layers, obj.ground_res, nrXY, true, obj.int_th, obj.sf.localNodes);
        
            % pillar then taper upwards
            obj.mesh.addLayeredQuarterCylinder([0, 0, 0], Rbase, obj.pillar_layers, obj.pillar_res, nrXY, true, obj.int_th, obj.sf.localNodes);
            obj.mesh.nodes = obj.mesh.coneTransformationX(pillar_height, Rbase, Rtop, obj.mesh.nodes);

           obj.mesh.addPipe3D([0,0,-obj.ground_layers(end)] ,Rbase,2*Rbase,0,90,-obj.ground_layers(end),0,10,12,obj.ground_res(end),obj.sf.localNodes);
        end

        
    end
end