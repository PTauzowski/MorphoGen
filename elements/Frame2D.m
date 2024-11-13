classdef Frame2D < FiniteElement
    properties
        E,A,J,m
    end
      
    methods
        function obj = Frame2D(elems)
            obj = obj@FiniteElement(ShapeFunctionsFrame2D, elems);
             obj.eDofs=["ux" "uy" "fiz"];
             obj.results.names = ["N" "Ty" "Mz" ];
             obj.results.descriptions = ["normal force" "shear force" "bending moment" ];
                 
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

        function K = computeLocalStifnessMatrix(obj, nodes, varargin)
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
                  K(:,:,k) =  L(:,:,k)' * [EA/l   , 0.0      , 0.0      , -EA/l  , 0.0       , 0.0;...
		                0.0    , 12*EJ/l3 , -6*EJ/l2 , 0.0    , -12*EJ/l3 , -6*EJ/l2;...
		                0.0    , -6*EJ/l2 , 4*EJ/l   , 0.0    , 6*EJ/l2   , 2*EJ/l;...
		                -EA/l  , 0.0      , 0.0      , EA/l   , 0.0       , 0.0 ;...
		                0.0    , -12*EJ/l3, 6*EJ/l2  , 0.0    , 12*EJ/l3  , 6*EJ/l2;...
		                0.0    , -6*EJ/l2 , 2*EJ/l   , 0.0    , 6*EJ/l2   , 4*EJ/l   ]  * L(:,:,k);
             end
        end

        function K = computeStifnessMatrix(obj, nodes, varargin)
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              K  = zeros( dim , dim, nelems );
              Kl = obj.computeLocalStifnessMatrix(nodes,varargin);
              L  = obj.computeTransformationMatrix(nodes);              
              for k=1:nelems
                 K(:,:,k) = L(:,:,k)' * Ke * L(:,:,k);
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
                 K(:,:,k) =  [ 0  0  0  0  0  0;...
		                       0  0  0  0  0  0; ...
		                       0  0  0  0  0  0;...
		                       0  0  0  0  0  0;...
		                       0  0  0  0  0  0;...
		                       0  0  0  0  0  0 ];
             end
             K=K(:);
        end

        function plot(obj, nodes)
            hold on;
            daspect([1 1 1]);
            plot([ nodes(obj.elems(:,1),1) nodes(obj.elems(:,2),1) NaN(21,1) ]',...
                 [ nodes(obj.elems(:,1),2) nodes(obj.elems(:,2),2) NaN(21,1) ]',...
                    "LineStyle","-","Marker","o","Color","k","LineWidth",2);
        end
        function plotMap(nodes,q,mapName,scd)
        end

        function N = shapeMatrix( obj, points )
        end
        function P = loadLineIntegral(obj, mode, nodes, fedges, dofnames, valueFn )
        end

        function initializeResults(obj)
        end

        function computeResults(obj,nodes,q,varargin)
        end

        function plotWired(obj,nodes,varargin)
            obj.plot(nodes);
        end


    end
end

