classdef DOFManagerUniform < DOFManager
    
    methods
        function  [I,J,V,Ksize] = getIndices( obj, fElem )
            nelems = size( fElem.elems, 1 );
            nnodes = size( fElem.elems, 2 );
            ndim=size(obj.eDofs,2);
            Kdim  = ndim * nnodes;
            Ksize = Kdim * Kdim;
            [ix, iy] = meshgrid( 1:Kdim, 1:Kdim );
            alldofs = (repelem( obj.elems, 1, ndim)-1)*ndim+repmat((1:ndim)',nelems,nnodes);
            I = alldofs(1:nelems,ix(:));
            J = alldofs(1:nelems,iy(:));
            V = alldofs;
        end
    end
end

