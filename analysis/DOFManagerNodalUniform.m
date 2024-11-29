classdef DOFManagerNodalUniform < DOFManager 

    properties
        nDOFs;
    end

    methods   
        function obj = DOFManagerNodalUniform(fElems)
            totalDofs=fElems{1}.eDofs;
            for k=2:size(fElems,2)
                totalDofs=[totalDofs fElems{k}.eDofs];               
            end
            obj.nDOFs = unique(totalDofs, 'stable');
        end

        function  [I,J,V,Ksize] = getIndices( obj, fElem )
            nelems = size( obj.elems, 1 );
            nnodes = size( obj.elems, 2 );
            [~,~,idofs] = intersect(fElem.eDofs, obj.nDOFs);
            Kdim  = size(fElem.eDofs,2) * nnodes;
            Ksize = Kdim * Kdim;
            [ix, iy] = meshgrid( 1:Kdim, 1:Kdim );
            alldofs = (repelem( obj.elems, 1, size(fElem.eDofs,2))-1)*size(obj.nDOFs,2)+repmat(idofs',nelems,nnodes);
            I = alldofs(1:nelems,ix(:));
            J = alldofs(1:nelems,iy(:));
            V = alldofs;
        end

    end
end

