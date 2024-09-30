classdef FEModel < handle
    properties
            felems, nodes, ndofs, supports, P, q, rotations;
            selTolerance;
    end

    methods
        function obj = FEModel(felems, nodes)
            if iscell(felems)
                obj.felems = felems;
            else
                obj.felems = { felems };
            end
            obj.nodes = nodes;
            obj.ndofs = obj.felems{1}.ndofs;
            for k=max(size(obj.felems))
                obj.ndofs = union(obj.ndofs,obj.felems{k}.ndofs);
            end
            obj.selTolerance=1.0E-05;
            obj.Pnodal = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
            obj.Pfem = []; %zeros( size(mesh.nodes,1) * size(obj.ndofs,2), 1 );
            obj.qnodal = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
            obj.supports = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
        end
        
    end

end

