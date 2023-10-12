classdef TopologyOptimization99 < TopologyOptimization
    %TOPOLOGYOPTIMIZATION99 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nelx, nely, volfrac, penal, rmin, ke, I, J, K, F, U, fixeddofs, alldofs, freedofs,supports, eqs;
    end
    
    methods
        function obj = TopologyOptimization99(nelx,nely,volfrac,penal,rmin)
           obj=obj@TopologyOptimization(rmin);
           obj.nelx=nelx;
           obj.nely=nely;
           obj.volfrac=volfrac;
           obj.penal=penal;
           obj.rmin=rmin;
           obj.x(1:nely,1:nelx) = volfrac; 
           E = 1.; 
           nu = 0.3;
           k=[ 1/2-nu/6  1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8  nu/6  1/8-3*nu/8];
           obj.ke = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                                 k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                                 k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                                 k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                                 k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                                 k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                                 k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                                 k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
           obj.I=zeros(nelx*nely*64,1);
           obj.J=zeros(nelx*nely*64,1);
           obj.K=zeros(nelx*nely*64,1);
           ielem=1;
           for elx = 1:nelx
                for ely = 1:nely
                    n1 = (nely+1)*(elx-1)+ely; 
                    n2 = (nely+1)* elx   +ely;
                    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
                    Ie=repmat(edof,1,8)';
                    obj.I(ielem:ielem+64-1)=Ie(:);
                    obj.J(ielem:ielem+64-1)=repmat(edof,8,1);
                    ielem=ielem+64;
                end
           end
           obj.F = zeros(2*(obj.nely+1)*(obj.nelx+1),1); 
           obj.U = zeros(2*(obj.nely+1)*(obj.nelx+1),1);
           obj.supports = zeros(2*(obj.nely+1)*(obj.nelx+1),1);
           obj.F(2,1) = -1;
           obj.fixeddofs   = union([1:2:2*(obj.nely+1)],[2*(obj.nelx+1)*(obj.nely+1)]);
           obj.alldofs     = [1:2*(obj.nely+1)*(obj.nelx+1)];
           obj.freedofs    = setdiff(obj.alldofs,obj.fixeddofs);
           obj.supports(obj.fixeddofs,:)=1;
           obj.eqs=LinearEquationsSystem(obj.I,obj.J,obj.supports);
           obj.createFilteringProperties3()
        end
        
        function x = solve(obj)
          % INITIALIZE
            x(1:obj.nely,1:obj.nelx) = obj.volfrac; 
            loop = 0; 
            change = 1.;
    
% START ITERATION
            while change > 0.003 
                loop = loop + 1;
                xold = x;
% FE-ANALYSIS
                [U]=obj.fe(x);         
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
                c = 0.;
                for ely = 1:obj.nely
                  for elx = 1:obj.nelx
                    n1 = (obj.nely+1)*(elx-1)+ely; 
                    n2 = (obj.nely+1)* elx   +ely;
                    Ue = obj.U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
                    c = c + x(ely,elx)^obj.penal*Ue'*obj.ke*Ue;
                    dc(ely,elx) = -obj.penal*x(ely,elx)^(obj.penal-1)*Ue'*obj.ke*Ue;
                  end
              end
% FILTERING OF SENSITIVITIES
              %[dc]   = obj.check(obj.nelx,obj.nely,obj.rmin,x,dc); 
              dc(:)   = obj.check3(dc(:));
              x(:)   = obj.check3(x(:));
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
              [x,lmid]    = obj.oc(obj.nelx,obj.nely,x,obj.volfrac,dc); 
% PRINT RESULTS
              change = max(max(abs(x-xold)));
              disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
                    ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(obj.nelx*obj.nely)) ...
                    ' ch.: ' sprintf('%6.3f',change ) ...
                    ' lmid.: ' sprintf('%6.3f',lmid )])
              colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
% PLOT DENSITIES  
            end  

            x=obj.x;
        end

        function [U]=fe(obj,x)
       
            ielem=1;
            for elx = 1:obj.nelx
                for ely = 1:obj.nely
                    obj.K(ielem:ielem+64-1) = x(ely,elx)^obj.penal*obj.ke;
                    ielem=ielem+64;
                end
            end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
            
% SOLVING
            obj.U(obj.freedofs,:) = obj.eqs.solveSimple(obj.K,obj.F);      
            obj.U(obj.fixeddofs,:)= 0;
            U=obj.U;
        end

        function [dcn]=check(obj,nelx,nely,rmin,x,dc)
            dcn=zeros(nely,nelx);
            for i = 1:nelx
                for j = 1:nely
                    sum=0.0; 
                    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
                        for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                            fac = rmin-sqrt((i-k)^2+(j-l)^2);
                            sum = sum+max(0,fac);
                            dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
                        end
                    end
                    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
                end
            end
        end
        
        function [xnew,lmid]=oc(obj,nelx,nely,x,volfrac,dc)  
            l1 = 0; l2 = 10000; move = 0.2;
            while (l2-l1 > 1e-4)
              lmid = 0.5*(l2+l1);
              xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
              if sum(sum(xnew)) - volfrac*nelx*nely > 0
                l1 = lmid;
              else
                l2 = lmid;
              end
            end

        end
        function createFilteringProperties4(obj)
            tne = obj.nelx*obj.nely;
            celemsi = cell( tne, 2 );
            xs = zeros( tne, 2 );
            k=1;
            for i=1:obj.nelx
                for j=1:obj.nely
                    xs(k, :) = [j+0.5 i+0.5];
                    k=k+1;
                end
            end
            np = size(xs,1);
            [Xsort, Ix] = sort(xs(:,1));
            [Ysort, Iy] = sort(xs(:,2));
            rIx(Ix) = 1:np;
            rIy(Iy) = 1:np;
            Ivx = obj.intervals( Xsort );
            Ivy = obj.intervals( Ysort );
            k=1;
            rf=obj.Rf;
            for k=1:size(xs,1)
                ip=intersect( Ix(Ivx(rIx(k),1):Ivx(rIx(k),2)), Iy(Ivy(rIy(k),1):Ivy(rIy(k),2)) );
                [~,d] = dsearchn( xs(k,:), xs(ip,:));
                celemsi{ k, 1 } = ip( d <= rf );
                celemsi{ k, 2 } = 1 - d( d <= rf ) / rf;
                k=k+1;
            end
            obj.celems=celemsi;
        end
        function elemValue = check4(obj,elemValue)
            dcn    = elemValue;
            dcn(:) = 0;
            for k=1:size(elemValue,1)
                fac = sum( obj.celems{k,2}(:) );
                if ( size( obj.celems{k}, 1) == 0 )
                    dcn( k ) = elemValue(k);
                else
                    dcn( k ) = sum(obj.celems{k, 2}(:).*elemValue(obj.celems{k, 1}(:)))/fac;
                end
            end
            elemValue = dcn;
        end
        function createFilteringProperties3(obj)
            tne = obj.nelx*obj.nely;
            centroids = zeros( tne, 2 );
            k=1;
            for i=1:obj.nelx
                for j=1:obj.nely
                    centroids(k, :) = [j+0.5 i+0.5];
                    k=k+1;
                end
            end

            weights   = sparse(tne, tne);
            distances = zeros( tne, 1);
            Rf = obj.Rf;
            for i=1:size( centroids, 1)
                    for j=1:size( centroids, 1)
                        distances(j) = norm( centroids(i, :) - centroids(j, :));
                    end
                    neighbours = find( distances <= Rf ); 
                    sum_near_dist = 0;
                    for k=1:max(size(neighbours))
                        sum_near_dist = sum_near_dist + (Rf-distances( neighbours(k)));
                    end
                    for j=1:max(size(neighbours))
                        weights(i, neighbours(j)) = (Rf-distances( neighbours(j)))/sum_near_dist;
                    end
            end
            obj.distances = distances;
            obj.weights = weights;
        end
        
        function elemValue = check3(obj,elemValue)
            elemValue = obj.weights * elemValue;
        end

        function resp = finishCondition(x)
            resp = sum(x);
        end

    end
end

