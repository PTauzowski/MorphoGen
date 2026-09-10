classdef FEModel < handle
    properties
            fElems, mesh, dofTypes, modelDofs, dofToNodes, nodesToDofs, supports, P, q, rotations;
            selTolerance;
    end

    properties (Access=protected)
            % Density/design field. Carried here to match the pre-refactor
            % FEModel, which ModelLinear and every model fixture rely on.
            % Conceptually it belongs on ModelLinear, its only consumer —
            % relocating it is a change the merge does not force (Rule 2).
            x;
    end

    methods

        function obj = FEModel( felems, mesh )
            % Model subclasses (ModelLinear and its fixtures) are built empty
            % and populate mesh/fElems in their own constructors, so they reach
            % this with no arguments. Only initialise DOFs when given both.
            if nargin < 2
                return;
            end
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
            % ============================================================
            % initDOFs  Initialize global DOF numbering and mapping tables
            %
            % This function constructs the mapping between mesh nodes and
            % degrees of freedom (DOFs) for all finite element (FE) classes.
            %
            % It supports multiple FE types coexisting on the same mesh
            % (e.g., solids, beams, shells), where nodes may carry different
            % sets of DOFs depending on which element classes touch them.
            %
            % Results stored in:
            %   obj.modelDofs   – list of all global DOF labels (e.g., {'ux','uy',...})
            %   obj.dofToNodes  – node number associated with each DOF
            %   obj.dofTypes    – unique list of all DOF kinds used in model
            %   obj.nodesToDofs – mapping: node → list of global DOF indices
            % ============================================================
        
            % Boolean matrix marking which nodes belong to each FE class
            elemsToNodes = false(obj.mesh.getNumberOfNodes(), obj.getFEClassesNumber());
        
            % Cell arrays for collecting DOF names and node numbers per node
            nodalDofs = repmat({''}, obj.mesh.getNumberOfNodes(), 1);
            nodeNums  = cell(obj.mesh.getNumberOfNodes(), 1);
        
            % List to collect all DOF types encountered in the model
            totalDofs = [];
        
            % ------------------------------------------------------------
            % Loop over all finite element classes in the model
            % ------------------------------------------------------------
            for k = 1:obj.getFEClassesNumber()
        
                % Retrieve node indices and DOF names for current FE class
                % Example: nodes = [1 2 3 4]; nDofs = {{'ux','uy','uz'}, ...}
                [nodes, nDofs] = obj.fElems{k}.getNodesDOFs();
        
                % Mark which nodes are used by this FE class
                elemsToNodes(obj.fElems{k}.elems, k) = true;
        
                % Collect all DOF kinds for global type registry
                totalDofs = [totalDofs obj.fElems{k}.eDofs];
        
                % --------------------------------------------------------
                % Merge DOFs contributed by this FE class into node lists
                % --------------------------------------------------------
                for i = 1:numel(nodes)
                    % Union ensures DOF kinds are unique per node
                    nodalDofs{nodes(i)} = union(nodalDofs{nodes(i)}, nDofs{i}, 'stable');
        
                    % Store node numbers (repeated for each DOF kind)
                    nodeNums{nodes(i)} = repelem(nodes(i), 1, numel(nodalDofs{nodes(i)}));
                end
            end
        
            % ------------------------------------------------------------
            % Remove dummy first entries ('') created by initialization
            % ------------------------------------------------------------
            uDofs = cellfun(@(a) a(2:end), nodalDofs, 'UniformOutput', false);
            nDofs = cellfun(@(a) a(2:end), nodeNums,  'UniformOutput', false);
        
            % Flatten nested cell arrays into long column vectors
            obj.modelDofs  = [uDofs{:}]';  % All DOF names, e.g. {'ux'; 'uy'; 'uz'; ...}
            obj.dofToNodes = [nDofs{:}]';  % Corresponding node numbers
        
            % Global registry of DOF types (e.g. ux, uy, uz, rx, ry, rz)
            obj.dofTypes = unique(totalDofs, 'stable');
        
            % Alias for clarity in next stage
            nodesToDofs = nDofs;
        
            % ------------------------------------------------------------
            % Assign contiguous global DOF indices node-by-node
            % ------------------------------------------------------------
            i = 1;
            for k = 1:obj.mesh.getNumberOfNodes()
                lDofs = numel(nodesToDofs{k});     % Number of DOFs on this node
                nodesToDofs{k} = i:i+lDofs-1;      % Assign global indices
                i = i + lDofs;                     % Increment global counter
            end
        
            % Store final node → DOF index map
            obj.nodesToDofs = nodesToDofs;
        end


        % function initDOFs(obj)
        %     elemsToNodes=false( obj.mesh.getNumberOfNodes(), obj.getFEClassesNumber() );
        %     nodalDofs=repmat({''},obj.mesh.getNumberOfNodes(),1);
        %     nodeNums=cell(obj.mesh.getNumberOfNodes(),1);
        %     totalDofs=[];
        %     for k=1:obj.getFEClassesNumber()
        %         [nodes nDofs] = obj.fElems{k}.getNodesDOFs();
        %         elemsToNodes(obj.fElems{k}.elems,k)=true;
        %         totalDofs=[totalDofs obj.fElems{k}.eDofs];
        %         for i = 1:size(nodalDofs(nodes),1)
        %             nodalDofs{nodes(i)} = union(nodalDofs{nodes(i)}, nDofs{i}, 'stable');
        %             nodeNums{nodes(i)} = repelem( nodes(i), 1,  numel(nodalDofs{nodes(i)}) );
        %         end
        %     end
        %     uDofs = cellfun(@(a) a(2:end), nodalDofs, 'UniformOutput', false);
        %     nDofs = cellfun(@(a) a(2:end), nodeNums, 'UniformOutput', false);
        %     obj.modelDofs = [uDofs{:}]';
        %     obj.dofToNodes = [nDofs{:}]';
        %     obj.dofTypes = unique(totalDofs, 'stable');
        %     nodesToDofs = nDofs;
        %     i=1;
        %     for k=1:obj.mesh.getNumberOfNodes()
        %         lDofs=numel(nodesToDofs{k});
        %         nodesToDofs{k}=i:i+lDofs-1;
        %         i=i+lDofs;
        %     end
        %     obj.nodesToDofs=nodesToDofs;
        % end

        function allocVectors = getAllocationVectors( obj, inds )
            a=obj.nodesToDofs(inds);
            allocVectors = obj.nodesToDofs(inds);
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

