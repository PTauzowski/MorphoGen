classdef FEModel < handle
    properties
            fElems, mesh, dofTypes, modelDofs, dofToNodes, nodesToDofs, supports, P, q, rotations;
            selTolerance;
    end

    methods

        function obj = FEModel( felems, mesh )
            obj.fElems = felems;
            obj.mesh = mesh;
            obj.initDOFs();
        end

        function ne = getFEClassesNumber(obj)
            ne=size(obj.fElems,2);
        end

        function ten = getFEInstancesNumber(obj)
            ne=[];
            ten=sum(cellfun( @(fe) [ne size(fe.elems,1)],obj.fElems));       
        end

        function initDOFs(obj)
            elemsToNodes=false( obj.mesh.getNumberOfNodes(), obj.getFEClassesNumber() );
            nodalDofs=repmat({''},obj.mesh.getNumberOfNodes(),1);
            nodeNums=cell(obj.mesh.getNumberOfNodes(),1);
            totalDofs=[];
            for k=1:obj.getFEClassesNumber()
                [nodes nDofs] = obj.fElems{k}.getNodesDOFs();
                elemsToNodes(obj.fElems{k}.elems,k)=true;
                totalDofs=[totalDofs obj.fElems{k}.eDofs];
                for i = 1:size(nodalDofs(nodes),1)
                    nodalDofs{nodes(i)} = union(nodalDofs{nodes(i)}, nDofs{i}, 'stable');
                    nodeNums{nodes(i)} = repelem( nodes(i), 1,  numel(nodalDofs{nodes(i)}) );
                end
            end
            uDofs = cellfun(@(a) a(2:end), nodalDofs, 'UniformOutput', false);
            nDofs = cellfun(@(a) a(2:end), nodeNums, 'UniformOutput', false);
            obj.modelDofs = [uDofs{:}]';
            obj.dofToNodes = [nDofs{:}]';
            obj.dofTypes = unique(totalDofs, 'stable');
            nodesToDofs = nDofs;
            i=1;
            for k=1:obj.mesh.getNumberOfNodes()
                lDofs=numel(nodesToDofs{k});
                nodesToDofs{k}=i:i+lDofs-1;
                i=i+lDofs;
            end
            obj.nodesToDofs=nodesToDofs;
        end

        function allocVectors = getAllocationVectors( obj, inds )
            a=obj.nodesToDofs(inds);
            allocVectors = obj.nodesToDofs(inds);
        end

        function fixDOF( obj, nodeSelector, dofs,  values )
            nodesToFix = nodeSelector.select(obj.mesh.nodes)
            
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

        function draw(obj, varargin)
             % Parse varargin
             daspect([1 1 1]);
             nodeMarker = '.';
             nodeColor  = 'b';
             elemStyle  = 'solid';
             elemColor  = [0.8 0.8 0.8];
             edgeStyle  = 'solid';
             edgeColor  = 'k';
             mode       = 'faceContour';
             if (nargin>1)
                 args=varargin;
                 for k = 1:2:length(varargin)
                    name = args{k};
                    value = args{k+1};
                    switch name
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
                    nElems = size(obj.fElems{k}.elems,1);
                    pedges=reshape(obj.fElems{k}.elems(:,obj.fElems{k}.shapeFn.edges)',edgesPerElem, nElems*size(obj.fElems{k}.shapeFn.edges,1))';
                    %pedgesSorted= sort(reshape(pedges',size(pedges,2)/edgesPerElem,size(pedges,1)*edgesPerElem)',2);
                    pedgesSorted= sort(pedges,2);
                    [uniqueColumns, ia, ic] = unique(pedgesSorted, 'rows', 'stable');
                    pedges=pedges(ia,:);
                    X=[reshape(obj.mesh.nodes(pedges,1),size(pedges)) repelem(NaN,size(ia,1),1)]';
                    Y=[reshape(obj.mesh.nodes(pedges,2),size(pedges)) repelem(NaN,size(ia,1),1)]';
                    line(X(:),Y(:),'Color' ,edgeColor);
                 end
                 if isprop(obj.fElems{k}.shapeFn,'fcontours')
                 end
             end
             switch obj.mesh.getDim()
                case 2
                     line(obj.mesh.nodes(:,1),obj.mesh.nodes(:,2),"LineStyle","none","Marker",nodeMarker,'MarkerFaceColor',nodeColor,'Color' ,nodeColor);
                case 3
                     line(obj.mesh.nodes(:,1),obj.mesh.nodes(:,2),obj.mesh.nodes(:,3),"LineStyle","none","Marker",nodeMarker,'MarkerFaceColor',nodeColor,'Color' ,nodeColor);                 
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

