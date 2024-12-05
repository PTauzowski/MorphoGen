classdef Frame3D < FiniteElement
      
    properties
        E,A,G,Jy,Jz,Ks;
       
    end
      
    methods
        function obj = Frame3D(elems, E,A,G,Jy,Jz,Ks)
             obj = obj@FiniteElement(ShapeFunctionsFrame3D,elems);

             obj.ndofs=["ux" "uy" "uz" "fix" "fiy" "fiz"];
             obj.results.names = ["N" "Ty" "Tz" "Ms" "My" "Mz"];
             obj.results.descriptions = ["normal force" "shear force Ty" "shear force Tz" ...
                 "torsion moment" "bending moment My" "bending moment Mz" ];
              obj.E=E;
              obj.A=A;
              obj.G=G;
              obj.Jy=Jy;
              obj.Jz=Jz;
              obj.Ks=Ks;
                 
        end

        function L = computeTransformationMatrix(obj, nodes)
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              L = zeros( dim , dim, nelems );
              c0 =zeros(3,3);
              for k=1:nelems
                 l=norm(nodes(obj.elems(k,2),:)-nodes(obj.elems(k,1),:));
                 dl=nodes(obj.elems(k,2),:)-nodes(obj.elems(k,1),:);
                 dx=dl(1); 
                 dy=dl(2);
                 dz=dl(3);
                 l1=norm(dl(1:2));
                 if l1>1.0E-8
                        cj=[  dx/l         dy/l        dz/l; ...
                             -dy/l1        dx/l1       0; ... 
                             -dx/l*dz/l1  -dz/l*dy/l1  l1/l ];
                 else
                     cj = [ 0    0  dz/l; ...
                            0    1  0;
                           -dz/l 0, 0];
                 end
                 
                 L(:,:,k) =  [ cj c0 c0 c0; ...
                               c0 cj c0 c0; ...
                               c0 c0 cj c0; ...
                               c0 c0 c0 cj];
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
                 K(:,:,k) = L(:,:,k)' * Kl(:,:,k) * L(:,:,k);
              end
              K=K(:);
        end

        function K = computeLocalStifnessMatrix(obj, nodes, varargin)
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              K  = zeros( dim , dim, nelems );
              EA=obj.E*obj.A; EJy=obj.E*obj.Jy; EJz=obj.E*obj.Jz; GKs=obj.G*obj.Ks;
              for k=1:nelems
                 l=norm(nodes(obj.elems(k,2),:)-nodes(obj.elems(k,1),:));
                 l2=l^2;
                 l3=l^3;
               
		         Ke(1,1)   = EA/l;       
                 Ke(2,2)   = 12*EJz/l3; 
                 Ke(3,3)   = 12*EJy/l3; 
                 Ke(4,4)   = GKs/l;
		         Ke(5,5)   = 4*EJy/l;    
                 Ke(6,6)   = 4*EJz/l;   
                 Ke(7,7)   = EA/l;      
                 Ke(8,8)   = 12*EJz/l3; 
                 Ke(9,9)   = 12*EJy/l3;
		         Ke(10,10) = GKs/l;      
                 Ke(11,11) = 4*EJy/l; 
                 Ke(12,12) = 4*EJz/l;

		         Ke(5,3)  = -6*EJy/l2;  
                 Ke(6,2)  =  6*EJz/l2;  
                 Ke(7,1)  = -EA/l;
		         Ke(8,2)  = -12*EJz/l3; 
                 Ke(8,6)  = -6*EJz/l2;
		         Ke(9,3)  = -12*EJy/l3; 
                 Ke(9,5)  = 6*EJy/l2;
		         Ke(10,4) = -GKs/l;
		         Ke(11,3) = -6*EJy/l2; 
                 Ke(11,5) = 2*EJy/l;  
                 Ke(11,9) = 6*EJy/l2;
		         Ke(12,2) =  6*EJz/l2; 
                 Ke(12,6) = 2*EJz/l;  
                 Ke(12,8) = -6*EJz/l2;

                 Ke(3,5)  = -6*EJy/l2;  
                 Ke(2,6)  =  6*EJz/l2;  
                 Ke(1,7)  = -EA/l;
		         Ke(2,8)  = -12*EJz/l3; 
                 Ke(6,8)  = -6*EJz/l2;
		         Ke(3,9)  = -12*EJy/l3; 
                 Ke(5,9)  = 6*EJy/l2;
		         Ke(4,10) = -GKs/l;
		         Ke(3,11) = -6*EJy/l2; 
                 Ke(5,11) = 2*EJy/l;  
                 Ke(9,11) = 6*EJy/l2;
		         Ke(2,12) =  6*EJz/l2; 
                 Ke(6,12) = 2*EJz/l;  
                 Ke(8,12) = -6*EJz/l2;
                 K(:,:,k) = Ke;
              end
        end

        function plotMap(obj)
        end
        function plot (obj, nodes)    
            hold on;
            daspect([1 1 1]);
            nelems=size(obj.elems,1);
             plot3([ nodes(obj.elems(:,1),1) nodes(obj.elems(:,2),1) NaN(nelems,1) ]',...
                   [ nodes(obj.elems(:,1),2) nodes(obj.elems(:,2),2) NaN(nelems,1) ]',...
                   [ nodes(obj.elems(:,1),3) nodes(obj.elems(:,2),3) NaN(nelems,1) ]',...
                    "LineStyle","-","Marker","o","Color","k","LineWidth",2);
        end
        function plotWired (obj)
        end
        function [Fel, Feg] = computeResults(obj,nodes,qnodal) 
            L  = obj.computeTransformationMatrix(nodes); 
            Ke = obj.computeLocalStifnessMatrix(nodes);
            Fel = zeros(12,size(L,3));
            Feg = zeros(12,size(L,3));
            nnodes=2;
            ndofs=6;
            nelems=size(obj.elems,1);
            qelems = reshape( qnodal( obj.elems',:)', nnodes * ndofs, nelems );
            for k=1:size(L,3)
                Fel(:,k) = Ke(:,:,k)*L(:,:,k)*qelems(:,k);
                Feg(:,k) = L(:,:,k)'*Fel(:,k);
            end
        end
        function initializeResults(obj)	
        end
        function loadLineIntegral(obj)
        end
            
        function shapeMatrix (obj)     	
        end
        
        end
end


