classdef FEModel < handle
    properties
            fElems, nodes, modelDofs, supports, P, q, rotations;
            selTolerance;
    end

    methods
    end

    methods (Access=private)

        function obj = setGeometry(fElems, nodes)
            if iscell(fElems)
                obj.fElems = fElems;
            else
                obj.fElems = { fElems };
            end
            obj.nodes = nodes;
            obj.modelDofs = obj.fElems{1}.elemDofs;
            for k=max(size(obj.fElems))
                obj.modelDofs = union(obj.modelDofs,obj.fElems{k}.eDofs);
            end
            obj.selTolerance=norm(max(nodes),2)*1.0E-06;
        end
        
        
    end

end

