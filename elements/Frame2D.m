classdef Frame2D < FiniteElement
    properties
        E,A,J,m
    end
      
    methods
        function obj = Frame2D(elems)
            obj = obj@FiniteElement(ShapeFunctionsFrame2D, elems);
        end
        function L = computeTransformationMatrix(obj, nodes)
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              L = zeros( dim , dim, nelems );
             for k=1:nelems
                 l=norm(nodes(obj.elems(k,2))-nodes(obj.elems(k,1)));
                 c=nodes(obj.elems(k,2),1)-nodes(obj.elems(k,1),1)/l;
                 s=nodes(obj.elems(k,2),2)-nodes(obj.elems(k,1),2)/l;
                 L(:,:,k) =  [ c, s, 0 , 0 ,0, 0;...
		                      -s, c, 0 , 0 ,0, 0; ...
		                       0, 0, 1 , 0, 0, 0;...
		                       0, 0, 0,  c, s, 0 ;...
		                       0, 0, 0, -s, c, 0;...
		                       0, 0, 0,  0, 0, 1 ];
             end
        end

        function K = computeStifnessMatrix(obj, nodes, varargin)
              EA=obj.E*obj.A;
              EJ=obj.E*obj.J;
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              K = zeros( dim , dim, nelems );
              L = obj.computeTransformationMatrix(nodes);
             for k=1:nelems
                 l=norm(nodes(obj.elems(k,2))-nodes(obj.elems(k,1)));
                 l2=l^2;
                 l3=l^3;
                  K(:,:,k) =  L(:,:,k)' * [EA/l   , 0.0      , 0.0      , -EA/l  , 0.0       , 0.0;...
		                0.0    , 12*EJ/l3 , -6*EJ/l2 , 0.0    , -12*EJ/l3 , -6*EJ/l2;...
		                0.0    , -6*EJ/l2 , 4*EJ/l   , 0.0    , 6*EJ/l2   , 2*EJ/l;...
		                -EA/l  , 0.0      , 0.0      , EA/l   , 0.0       , 0.0 ;...
		                0.0    , -12*EJ/l3, 6*EJ/l2  , 0.0    , 12*EJ/l3  , 6*EJ/l2;...
		                0.0    , -6*EJ/l2 , 2*EJ/l   , 0.0    , 6*EJ/l2   , 4*EJ/l   ]  * L(:,:,k);
             end
             K=K(:);
        end
        function K = computeGeometricStifnessMatrix(obj, nodes, varargin)
              EA=obj.E*obj.A;
              EJ=obj.E*obj.J;
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              K = zeros( dim , dim, nelems );
             for k=1:nelems
                 l=norm(nodes(obj.elems(k,2))-nodes(obj.elems(k,1)));
                 l2=l^2;
                 l3=l^3;
                 K(:,:,k) =  [ 0, 0, 0, 0 ,0, 0;...
		                       0, 0, 0, 0 ,0, 0; ...
		                       0, 0, 0, 0, 0, 0;...
		                       0, 0, 0, 0, 0, 0;...
		                       0, 0, 0, 0, 0, 0;...
		                       0, 0, 0, 0, 0, 0 ];
             end
             K=K(:);
        end
    end
end

