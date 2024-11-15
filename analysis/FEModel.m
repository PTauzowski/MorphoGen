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

        function plotNodes(obj,marker, color)
            switch obj.mesh.getDim()
                case 2
                     line(obj.mesh.nodes(:,1),obj.mesh.nodes(:,2),"LineStyle","none","Marker",marker,'MarkerFaceColor',color);
                case 3
                     line(obj.mesh.nodes(:,1),obj.mesh.nodes(:,2),obj.mesh.nodes(:,3),"LineStyle","none","Marker",marker,'MarkerFaceColor',color);                 
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

        function draw(obj, varanargin)
             % Parse varargin

             nodeMarker = '.';
             nodeColor  = 'k';
             elemStyle  = 'solid';
             elemColor  = [0.8 0.8 0.8];
             edgeStyle  = 'solid';
             edgeColor  = 'k';
             mode       = 'faceContour';
             if (nargin>1)
                 args=varanargin;
                 for k = 1:2:length(varargin)
                    name = args{k};
                    value = args{k+1};
                    switch lower(name)
                        case 'nodeMarker'
                            nodeMarker = value;
                        case 'nodeColor'
                            nodeColor = value;
                        case 'elemStyle'
                            elemStyle = value;
                        case 'elemColor'
                            elemColor = value;
                        case 'edgeStyle'
                            edgeStyle = value;  
                        case 'edgeColor'
                            edgeColor = value; 
                        case 'mode'
                            marker = value;       
                        otherwise
                            error('Unknown plot style parameter: %s', name);
                    end
                 end
             end
             for k=1:size(obj.fElems,2)
                 if isprop(obj.fElems{k}.shapeFn, 'edges')
                    edgesPerElem = size(obj.fElems{k}.shapeFn.edges,2);
                    pedges=obj.fElems{k}.elems(1:2,obj.fElems{k}.shapeFn.edges);
                    pedges= sort(reshape(pedges',size(pedges,2)/edgesPerElem,size(pedges,1)*edgesPerElem)',2);
                    [uniqueColumns, ia, ic] = unique(pedges, 'rows', 'stable');
                    duplicateColumns = setdiff(1:size(pedges, 1), ia);
                    a=[uniqueColumns repelem(NaN,size(uniqueColumns,1))']';
                    edgesp=a(:);
                    line([obj.mesh.nodes(edgesp,1)],obj.mesh.nodes(edgesp,2),obj.mesh.nodes(edgesp,1))
                 end
                 if isprop(obj.fElems{k}.shapeFn,'fcontours')
                 end
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

