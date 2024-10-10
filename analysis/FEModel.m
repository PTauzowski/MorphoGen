classdef FEModel < handle
    properties
            fElems, mesh, modelDofs, supports, P, q, rotations;
            selTolerance;
    end

    methods

        function obj = FEModel()
            obj.fElems = {};
        end

        function addFiniteElements( obj, fElems )
            obj.fElems = { obj.fElems fElems };
            % if iscell(fElems)
            %     obj.fElems = { obj.fElems fElems };
            % else
            %     obj.fElems = { obj.fElems fElems };
            % end
        end

    end

    methods (Access=private)

        function obj = setModelDOFs(obj)
            nnodes=size(obj.mesh.nodes,1);
            obj.modelDof=cell(nnodes,1);
            for k=max(size(obj.fElems))
                for l=1:size(obj.fElems{k}.elems,1)
                    obj.modelDofs{obj.fElems{k}.elems(l,1)} = union(obj.modelDofs{obj.fElems{k}.elems(l,1)},obj.fElems{k}.eDofs);
                end
            end
            obj.selTolerance=norm(max(nodes),2)*1.0E-06;
        end
        
        
    end

end

