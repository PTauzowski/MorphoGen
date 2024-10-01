classdef DOFManager < handle
   
    
    properties
        dofmapToFEM, nnodeFromFEM, ndofs, gdofs;
        nnodes;
    end
    
    methods

        function obj = setUniformDOFS( obj, nnodes, ndofs )
            obj.gdofs=repelem(ndofs,1,nnodes);
            obj.ndofs=ndofs;
            nelemdofs=size(ndofs,1);
            obj.dofmapToFEM=zeros(nelemdofs,nnodes);
            obj.dofmapToFEM(:)=1:(nelemdofs*nnodes);
        end

        function dim = getDimension(obj)
            dim=max(size(obj.dofmapToFEM(:)));
        end

    end

end

