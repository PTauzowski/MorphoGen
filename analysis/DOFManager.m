classdef DOFManager < handle
   
    
    properties
        dofmapToFEM, ndofs, gdofs;
    end
    
    methods

        function obj = setUniformDOFS( obj, nnodes, ndofs )
            obj.gdofs=repelem(ndofs,1,nnodes);cell
            obj.ndofs=ndofs;
            nelemdofs=size(ndofs,2);
            obj.dofmapToFEM=zeros(nelemdofs,nnodes);
            obj.dofmapToFEM(:)=1:(nelemdofs*nnodes);            
            [x,y]=meshgrid(1:nnodes,1:nelemdofs);
        end

        function dim = getDimension(obj)
            dim=max(size(obj.dofmapToFEM(:)));
        end


    end

end

