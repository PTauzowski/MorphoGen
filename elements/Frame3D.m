classdef Frame3D < FiniteElement
      
    properties
        E,A,G,Jy,Jz,Ks,betas;
       
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

        function L = computeTransformationMatrixAnsys(obj, nodes)
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              L = zeros( dim , dim, nelems );
              T = zeros(3,3);
              T0 = zeros(3,3);
              for k=1:nelems
                 l=norm(nodes(obj.elems(k,2),:)-nodes(obj.elems(k,1),:));
                 dl=nodes(obj.elems(k,2),:)-nodes(obj.elems(k,1),:);
                 dx=dl(1); 
                 dy=dl(2);
                 dz=dl(3);
                 lxy=norm(dl(1:2));
                 d=0.0001*l;
                 x1 = nodes(obj.elems(k,1),1);
                 x2 = nodes(obj.elems(k,2),1);
                 y1 = nodes(obj.elems(k,1),2);
                 y2 = nodes(obj.elems(k,2),2);
                 z1 = nodes(obj.elems(k,1),3);
                 z2 = nodes(obj.elems(k,2),3);
                 if lxy>d
                    s1= (y2-y1)/lxy;  
                    c1= (x2-x1)/lxy;
                    
                 else
                     s1=0.0;
                     c1=1.0;
                 end
                 s2=(z2-z1)/l;
                 c2=lxy/l;
                 s3=0; %sin(deg2rad(obj.betas(k)));
                 c3=1; %cos(deg2rad(obj.betas(k)));

                 T= [  c1*c2               s1*c2              s2; ...
                      (-c1*s2*s3-s1*c3)   (-s1*s2*s3+c1*c3)   s3*c2; ...
                      (-c1*s2*c3-s1*s3)   (-s1*s2*c3-c1*s3)   c3*c2 ];
                 
                 L(:,:,k) =  [ T T0 T0 T0; ...
                               T0 T T0 T0; ...
                               T0 T0 T T0; ...
                               T0 T0 T0 T ];
             end
        end

        


        function K = computeStifnessMatrix(obj, nodes, varargin)
              nelems = size(obj.elems,1);
              nnodes = size(obj.elems,2);
              ndofs = size( obj.ndofs,2);
              dim = nnodes * ndofs;
              K  = zeros( dim , dim, nelems );
              Kl = obj.computeLocalStifnessMatrix(nodes,varargin);
              L = obj.computeTransformationMatrixAnsys(nodes);
              %L  = obj.computeTransformationMatrixAnsys2(nodes,'gamma_deg',obj.betas);
              %L  = obj.computeTransformationMatrixRotated(nodes,'gamma_deg',obj.betas); 
              %L  = obj.computeTransformationMatrix(nodes);  
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

        function plotMap(obj,nodes)
        end
        
        function plot(obj, nodes)    
            nelems=size(obj.elems,1);
            plot3([ nodes(obj.elems(:,1),1) nodes(obj.elems(:,2),1) NaN(nelems,1) ]',...
                   [ nodes(obj.elems(:,1),2) nodes(obj.elems(:,2),2) NaN(nelems,1) ]',...
                   [ nodes(obj.elems(:,1),3) nodes(obj.elems(:,2),3) NaN(nelems,1) ]',...
                    "LineStyle","-","Marker","o","MarkerEdgeColor",'r',"MarkerFaceColor",'r',"Color","k","LineWidth",3);
        end
        function plotSelected(obj, nodes, idx)    
            nelems=size(obj.elems,1);
            plot3([ nodes(obj.elems(idx,1),1) nodes(obj.elems(idx,2),1) NaN(size(idx,2),1) ]',...
                   [ nodes(obj.elems(idx,1),2) nodes(obj.elems(idx,2),2) NaN(size(idx,2),1) ]',...
                   [ nodes(obj.elems(idx,1),3) nodes(obj.elems(idx,2),3) NaN(size(idx,2),1) ]',...
                    "LineStyle","-","Marker",".","Color","m","LineWidth",5);
        end

        function h = plotLocalCS(obj, nodes, varargin)
           % Draw one triad per element at its midpoint.
            % Options:
            %   'gamma_deg'   : nelemsx1 or scalar axial twist (deg)
            %   'y_ref'       : nelemsx3 desired local Y vectors (global) – overrides gamma
            %   'scale'       : arrow length (default 0.25 * mean element length)
            %   'offsetLocal' : [ox oy oz] shift of base point in LOCAL coords (default [0 0 0])
            %   'idx'         : which elements to draw (default all)
            %   'ax'          : target axes (default gca)
            
                p = inputParser;
                addParameter(p,'gamma_deg',[],@isnumeric);
                addParameter(p,'y_ref',[],@(x)isnumeric(x) && (isempty(x) || size(x,2)==3));
                addParameter(p,'scale',[],@isnumeric);
                addParameter(p,'offsetLocal',[0 0 0],@(x)isnumeric(x)&&numel(x)==3);
                addParameter(p,'idx',[],@isnumeric);
                addParameter(p,'ax',[],@(x) isempty(x) || isgraphics(x,'axes'));
                parse(p,varargin{:});
                gamma_deg   = p.Results.gamma_deg;
                y_ref       = p.Results.y_ref;
                offsetLocal = p.Results.offsetLocal(:).';
                ax          = p.Results.ax; if isempty(ax), ax = gca; end
            
                ne = size(obj.elems,1);
                if isempty(p.Results.idx), idx = 1:ne; else, idx = p.Results.idx(:).'; end
                if isempty(gamma_deg), gamma_deg = zeros(ne,1); end
                if isrow(gamma_deg),   gamma_deg = gamma_deg(:); end
                if isscalar(gamma_deg), gamma_deg = repmat(gamma_deg, ne, 1); end
            
                lens = vecnorm(nodes(obj.elems(:,2),:) - nodes(obj.elems(:,1),:), 2, 2);
                sc   = p.Results.scale;
                if isempty(sc), sc = 0.25 * mean(lens); end
            
                nplot = numel(idx);
                P  = zeros(nplot,3); Ux = P; Uy = P; Uz = P;
            
                k = 0;
                for e = idx
                    i = obj.elems(e,1);  j = obj.elems(e,2);
                    p1 = nodes(i,:);     p2 = nodes(j,:);
                    mid = 0.5*(p1+p2);
            
                    ex = p2 - p1; ex = ex / norm(ex);
            
                    % local Y seed
                    if ~isempty(y_ref)
                        ey0 = y_ref(e,:);
                        ey0 = ey0 - dot(ey0,ex)*ex; ey0 = ey0 / norm(ey0);
                    else
                        if abs(ex(3)) < 0.95
                            ref = [0 0 1];
                        else
                            ref = [1 0 0];
                        end
                        ey0 = ref - dot(ref,ex)*ex; ey0 = ey0 / norm(ey0);
                    end
            
                    g  = deg2rad(gamma_deg(e));
                    ey =  ey0*cos(g) + cross(ex,ey0)*sin(g);
                    ez =  cross(ex,ey);
                    ey = ey / norm(ey); ez = ez / norm(ez);
            
                    base = mid + offsetLocal(1)*ex + offsetLocal(2)*ey + offsetLocal(3)*ez;
            
                    k = k+1;
                    P(k,:)  = base;
                    Ux(k,:) = sc * ex;
                    Uy(k,:) = sc * ey;
                    Uz(k,:) = sc * ez;
                end
            
                hold_state = ishold(ax); hold(ax,'on');
                hx = quiver3(ax, P(:,1),P(:,2),P(:,3), Ux(:,1),Ux(:,2),Ux(:,3), 0, 'Color',[0.85 0.1 0.1], 'LineWidth',1.2);
                hy = quiver3(ax, P(:,1),P(:,2),P(:,3), Uy(:,1),Uy(:,2),Uy(:,3), 0, 'Color',[0.00 0.6 0.0], 'LineWidth',1.2);
                hz = quiver3(ax, P(:,1),P(:,2),P(:,3), Uz(:,1),Uz(:,2),Uz(:,3), 0, 'Color',[0.10 0.1 0.9], 'LineWidth',1.2);
                if ~hold_state, hold(ax,'off'); end
                axis(ax,'equal');
            
                h = struct('x',hx,'y',hy,'z',hz);
         end


        function plotLocalCS_old(obj, nodes, scale, zoffset, relative)
            if nargin<3 || isempty(scale),   scale   = 0.15; end
            if nargin<4 || isempty(zoffset), zoffset = 0.0;  end
            if nargin<5 || isempty(relative),relative = false; end
        
            ne = size(obj.elems,1);
            C  = zeros(ne,3); EX=C; EY=C; EZ=C; Ls=zeros(ne,1);
        
            L  = obj.computeTransformationMatrixAnsys(nodes);
        
            for k=1:ne
                i = obj.elems(k,1); j = obj.elems(k,2);
                x1 = nodes(i,:); x2 = nodes(j,:);
                C(k,:)  = 0.5*(x1 + x2);
                Ls(k)   = norm(x2 - x1);
                R       = L(1:3,1:3,k);  % rows: ex; ey; ez
                EX(k,:) = R(1,:); EY(k,:) = R(2,:); EZ(k,:) = R(3,:);
            end
        
            % offset along LOCAL z if requested
            if relative, C = C + (zoffset.*Ls).*EZ; else, C = C + zoffset.*EZ; end
        
            holdState = ishold; hold on
            % x=red, y=green, z=blue
            quiver3(C(:,1),C(:,2),C(:,3), scale*EX(:,1), scale*EX(:,2), scale*EX(:,3), 0, 'Color',[0.85 0.10 0.10], 'LineWidth',1.3);
            quiver3(C(:,1),C(:,2),C(:,3), scale*EY(:,1), scale*EY(:,2), scale*EY(:,3), 0, 'Color',[0.10 0.60 0.10], 'LineWidth',1.3);
            quiver3(C(:,1),C(:,2),C(:,3), scale*EZ(:,1), scale*EZ(:,2), scale*EZ(:,3), 0, 'Color',[0.10 0.10 0.85], 'LineWidth',1.3);
            axis equal
            if ~holdState, hold off, end
        end
        
        
       
        function plotLocalForces(obj,nodes,qnodal,opts)
            if nargin<4, opts=struct; end
            opts = defaults(opts, 'scaleF',1e-3, 'scaleM',5e-4, 'sample','avg');
            
            [Fel,~] = obj.computeResults(nodes,qnodal);
            ne = size(obj.elems,1);
            
            C=zeros(ne,3); EX=C; EY=C; EZ=C;
            N=zeros(ne,1); Ty=N; Tz=N; Ms=N; My=N; Mz=N;
            
            L = obj.computeTransformationMatrix(nodes);
            for k=1:ne
                i=obj.elems(k,1); j=obj.elems(k,2);
                x1=nodes(i,:); x2=nodes(j,:); 
                C(k,:)  = 0.5*(x1+x2);
                R       = L(1:3,1:3,k);
                EX(k,:) = R(1,:); EY(k,:)=R(2,:); EZ(k,:)=R(3,:);
                f = Fel(:,k);
                switch opts.sample
                    case 'end1'
                        N(k)=f(1);  Ty(k)=f(2);  Tz(k)=f(3);
                        Ms(k)=f(4); My(k)=f(5);  Mz(k)=f(6);
                    otherwise % avg
                        N(k)=0.5*(f(1)-f(7));   Ty(k)=0.5*(f(2)-f(8));   Tz(k)=0.5*(f(3)-f(9));
                        Ms(k)=0.5*(f(4)-f(10)); My(k)=0.5*(f(5)-f(11));  Mz(k)=0.5*(f(6)-f(12));
                end
            end
            
            holdstate=ishold; hold on
            sF=opts.scaleF;
            quiver3(C(:,1),C(:,2),C(:,3), sF*N.*EX(:,1), sF*N.*EX(:,2), sF*N.*EX(:,3), 0, 'Color',[0.80 0 0],'LineWidth',1.6);
            quiver3(C(:,1),C(:,2),C(:,3), sF*Ty.*EY(:,1), sF*Ty.*EY(:,2), sF*Ty.*EY(:,3), 0, 'Color',[0 0.55 0],'LineWidth',1.6);
            quiver3(C(:,1),C(:,2),C(:,3), sF*Tz.*EZ(:,1), sF*Tz.*EZ(:,2), sF*Tz.*EZ(:,3), 0, 'Color',[0 0.2 0.8],'LineWidth',1.6);
            
            % moments as arcs -> one line per type, NaN-separated
            sM=abs(opts.scaleM);
            plot3(makeArcs(C,EX,sM*abs(Ms)),'-','Color',[0.85 0.45 0],'LineWidth',1.4);
            plot3(makeArcs(C,EY,sM*abs(My)),'-','Color',[0.45 0 0.85],'LineWidth',1.4);
            plot3(makeArcs(C,EZ,sM*abs(Mz)),'-','Color',[0 0.65 0.85],'LineWidth',1.4);
            
            axis equal
            if ~holdstate, hold off, end
            
                function [X,Y,Z]=makeArcs(C_,N_,R_)
                    nseg = 24; t = linspace(0,2*pi,nseg);
                    % basis in plane ⟂ normal
                    A = [1,0,0]; 
                    U = A - dot(A,N_,2).*N_; nz = vecnorm(U,2,2); 
                    swap = nz<1e-8; A(swap,:) = repmat([0,1,0],sum(swap),1);
                    U = A - dot(A,N_,2).*N_; U = U./vecnorm(U,2,2);
                    V = cross(N_,U,2);
                    % allocate and fill with NaNs between arcs
                    X = NaN(ne*nseg,1); Y=X; Z=X;
                    for p=1:ne
                        P = C_(p,:) + R_(p)*( U(p,:).*cos(t).' + V(p,:).*sin(t).' );
                        idx = (p-1)*nseg + (1:nseg);
                        X(idx)=P(:,1); Y(idx)=P(:,2); Z(idx)=P(:,3);
                    end
                end
        end


        function plotWired (obj)
        end

        function [Fel, Feg] = computeResults(obj,nodes,qnodal) 
            L = obj.computeTransformationMatrixAnsys(nodes);
            %L  = obj.computeTransformationMatrixAnsys2(nodes,'gamma_deg',obj.betas);
            %L  = obj.computeTransformationMatrixRotated(nodes,'gamma_deg',obj.betas);
            %L  = obj.computeTransformationMatrix(nodes); 
            Ke = obj.computeLocalStifnessMatrix(nodes);
            Fel = zeros(12,size(L,3));
            Feg = zeros(12,size(L,3));
            nnodes=2;
            ndofs=6;
            nelems=size(obj.elems,1);
            qelems = reshape( qnodal( obj.elems',:)', nnodes * ndofs, nelems );
            xend=nodes(end,:);
            Fels=zeros(6,size(nodes,1));
            for k=1:size(nodes,1)
                xk=nodes(k,:);
                Fels(1,k)=0;
                Fels(2,k)=0;
                Fels(3,k)=1;
                Fels(4,k)= (xend(2)-xk(2));
                Fels(5,k)=-(xend(1)-xk(1));
                Fels(6,k)=0;
            end
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


