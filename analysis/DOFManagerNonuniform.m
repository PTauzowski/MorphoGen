classdef DOFManagerNonuniform
  
    properties
        modelDofs, dofToNodes, dofTypes, nodesToDofs;
    end
    
    methods
        
        function obj = DOFManagerNonuniform(mesh,fElems)
            elemClassesNumber=size(fElems,2);
            elemsToNodes=false( mesh.getNumberOfNodes(), elemClassesNumber );
            nodeNums=cell(mesh.getNumberOfNodes(),1);
            totalDofs=[];
            for k=1:elemClassesNumber
                [nodes, nDofs] = fElems{k}.getNodesDOFs();
                elemsToNodes(fElems{k}.elems,k)=true;
            end
            [unique_connections, conn_classes, node_classes] = unique(elemsToNodes, 'rows');
            dofs_classes=cell(numel(conn_classes),1);
            dofs_numbers=zeros(numel(conn_classes),1);
            for k=1:numel(conn_classes)
                a=cellfun(@(a) [a.eDofs], {fElems{elemsToNodes(conn_classes(k),:)'}}, 'UniformOutput', false);
                dofs_classes{k}=unique([a{:}],'stable');
                dofs_numbers(k)=size(dofs_classes{k},2);
            end
            stringified = cellfun(@(x) strjoin(cellstr(x), ','), dofs_classes, 'UniformOutput', false);
            [~, dofs_classes_id, dofs_clasification] = unique(stringified);
            nodalDOFSClasses = dofs_clasification(node_classes);
            nodalDOFS = dofs_classes(nodalDOFSClasses);
            nodalDOFSnumber = dofs_numbers(nodalDOFSClasses);
            DOFsInds = [1; cumsum(nodalDOFSnumber)+1];
            globalDOFs = [ nodalDOFS{:} ]';
            for k=1:elemClassesNumber
                [nodes, nDofs] = fElems{k}.getNodesDOFs();
                elemDOFSclasses = nodalDOFSClasses( fElems{k}.elems );
                a = nodalDOFS(elemDOFSclasses);
                elemDOFs = [a{:}];
            end
        end

         function  [I,J,V,Ksize] = getIndices( obj, fElem )
            nelems = size( fElem.elems, 1 );
            nnodes = size( fElem.elems, 2 );
            a=obj.nodesToDofs(fElem.elems)';
            eDOFs=[a{:}];
            nDofs=obj.modelDofs(eDOFs);
            [newd,ai,idofs] = intersect(fElem.eDofs, nDofs,'stable');
            % for k=1:nelems
            %     a=obj.nodesToDofs(fElem.elems(1,:));
            %     eDofs=[a{:}];
            %     [~,~,idofs] = intersect(fElem.eDofs, obj.modelDofs(eDofs),'stable');
            % end
            Kdim  = size(eDOFs,2);
            Ksize = Kdim * Kdim;
            [ix, iy] = meshgrid( 1:Kdim, 1:Kdim );
            alldofs = eDOFs;
            I = alldofs(ix(:));
            J = alldofs(iy(:));
            V = alldofs;
        end
        
       
    end
end

