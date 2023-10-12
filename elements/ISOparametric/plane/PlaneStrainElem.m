classdef PlaneStrainElem < PlaneElem
    
    methods
        function obj = PlaneStrainElem(sf,p)
             obj = obj@PlaneElem(sf,p);
             obj.ndofs=["ux" "uy"];
        end
        
        function B = strainDerivMatrix( obj, dNx )
            ndofs = size(obj.elems,2) * size( obj.ndofs, 2 );
            nip   = size(dNx,3); 
            nnd = size(dNx,2); 
            %B = zeros( 3, ndofs, npg );
            B = cell(nip,1);
            for i=1:nip
                for j = 1:nnd
                      B{j}(1, 1:2:ndofs-1) = dNx(1,j,i);
                      B{j}(2, 2:2:ndofs)   = dNx(2,j,i);
                      B{j}(3, 1:2:ndofs-1) = dNx(2,j,i);
                      B{j}(3, 2:2:ndofs)   = dNx(1,j,i);
                end
            end
        end
        function setIsotropicMaterial( obj, E, nu, rho )
            D = E/(1+nu)/(1-2*nu)*[ 1-nu  nu 0; ...
                  nu 1-nu  0;
                  0  0  1-2*nu ];
            M =  diag([rho rho]); 
            obj.props.D = D;
            obj.props.M = M;
        end
    end
end

