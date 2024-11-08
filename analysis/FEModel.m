classdef FEModel < handle
    properties
            fElems, mesh, modelDofs, supports, P, q, rotations;
            selTolerance;
    end

    methods

        function obj = FEModel( felems, mesh )
            obj.fElems = felems;
            obj.mesh = mesh;
        end

        function addFiniteElements( obj, fElems )
            obj.fElems = { obj.fElems fElems };
        end

        function fixDOF( obj, nodeSelector, dofs,  values )
            nodesToFix = nodeSelector.select(obj.mesh.nodes);
            
        end

        function plotNodes(obj)
            switch obj.mesh.getDim()
                case 2
                     plot(obj.mesh.nodes(:,1),obj.mesh.nodes(:,2),"LineStyle","none","Marker","square",'MarkerFaceColor',[0.0,0.0,0.8]);
                case 3
                     plot(obj.mesh.nodes(:,1),obj.mesh.nodes(:,2),obj.mesh.nodes(:,3),"LineStyle","none","Marker","square",'MarkerFaceColor',[0.0,0.0,0.8]);                 
            end
        end

        function plotSelectedNodes(obj, selector ,col)
            selNodes = selector.select(obj.mesh.nodes);
            switch obj.mesh.getDim()
                case 2
                    plot(obj.mesh.nodes(selNodes,1),obj.mesh.nodes(selNodes,2),"LineStyle","none","Color",col,"Marker","square",'MarkerFaceColor',[0.5,0.0,0.5]);
                case 3
            end
        end

        function plot(obj)
            cellfun( @(fe) fe.plot(obj.mesh.nodes), obj.fElems);
        end

        function plotWired(obj)
            cellfun( @(fe) fe.plotWired(obj.mesh.nodes), obj.fElems);
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

