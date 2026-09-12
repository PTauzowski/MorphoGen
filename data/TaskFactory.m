classdef TaskFactory
    properties
        data
    end

    methods
        function obj = TaskFactory(data_filename)
            % If input is a filename ending with .json => load JSON
            if ischar(data_filename) || isstring(data_filename)
                [~, ~, ext] = fileparts(data_filename);
                if strcmpi(ext, '.json')
                    obj.data = jsondecode(fileread(data_filename));
                else
                    error('TaskFactory:InvalidFile', ...
                          'File must be a .json configuration file.');
                end
            elseif isstruct(data_filename)
                obj.data = data_filename;   % already-decoded JSON struct
            else
                error('TaskFactory:InvalidInput', ...
                      'Input must be a .json filename or a struct.');
            end
        end

        function mesh = createMesh(obj)
            switch obj.data.domain.shape
                case "rectangular"
                    l  = obj.data.domain.size.length;
                    h  = obj.data.domain.size.height;
                    sf = obj.getShapeFunction();
                    nx = obj.data.domain.mesh.nelx;
                    ny = obj.data.domain.mesh.nely;
                    mesh = Mesh();
                    mesh.addRectMesh2D(0, 0, l, h, nx, ny, sf.pattern);
            end
        end

        function sf = getShapeFunction(obj)
            switch obj.data.domain.shape_fn
                case "Q4"
                    sf = ShapeFunctionQ4;
                otherwise
                    error("Unknown shape function: %s", obj.data.shapeFunction.type);
            end
        end

        function fe = createFiniteElement(obj, mesh)
            fe = PlaneStressElem(obj.getShapeFunction(), mesh.elems);
        end

        function mat = createMaterial(obj)
            m = obj.data.material;
            mat = PlaneStressMaterial('mat1');
            mat.setElasticIzo(m.E, m.nu);
            mat.setMassIzoMatrix(m.rho);
        end

        function analysis = createAnalysis(obj, fe, mesh)
            analysis = LinearElasticityWeighted(fe, mesh, false);
            % ... configure supports etc ...
        end

        function topOpt = createTopOpt(obj)
            mesh = obj.createMesh();
            fe   = obj.createFiniteElement(mesh);
            fe.setMaterial(obj.createMaterial());
            analysis = obj.createAnalysis(fe, mesh);

            o = obj.data.topology;
            g = obj.data.geometry;

            Rfilter = obj.data.filter.RfilterFactor * g.h / g.res;

            topOpt = StressIntensityTopologyOptimizationVol( ...
                Rfilter, analysis, o.cutThreshold, o.penal, min(o.volumeFractions), true );
        end
    end
end
