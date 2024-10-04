classdef DOFManager < handle
   
    
    properties
        dofmapToFEM, ndofs, gdofs;
        nnodes;
    end
    
    methods

        function obj = setUniformDOFS( obj, nnodes, ndofs )
            obj.gdofs=repelem(ndofs,1,nnodes);
            obj.ndofs=ndofs;
            nelemdofs=size(ndofs,2);
            obj.dofmapToFEM=zeros(nelemdofs,nnodes);
            obj.dofmapToFEM(:)=1:(nelemdofs*nnodes);            
            [x,y]=meshgrid(1:nnodes,1:nelemdofs);
            obj.nnodeFromFEM=x(:);
            obj.dofFromFEM=y(:);
        end

        function dim = getDimension(obj)
            dim=max(size(obj.dofmapToFEM(:)));
        end


    end

end

