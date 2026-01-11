classdef Mesh < handle
    
    properties
        nodes, elems;
        tolerance = 10000;
        convex_hull_nodes;
    end
    
    methods
        function el2 = merge( obj, newNodes, newElems )
            [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
            newNodes = newNodes( si1, : );
            newElems = si2(newElems);
            if  size(newElems,2)==1
                newElems=newElems';
            end
           if size( obj.nodes, 1 ) == 0 
                  obj.nodes = newNodes;
                  obj.elems = newElems;
                  el2 = newElems;
           else
                [~,i1,i2] = intersect( round(obj.nodes .* obj.tolerance), round(newNodes.*obj.tolerance), 'rows' );
                nidx  = 1:size(newNodes,1);
                nidx( i2 ) = [];
                noi   = 1:size(nidx,2);
                noi = noi + size(obj.nodes,1);
                obj.nodes = [ obj.nodes; newNodes(nidx,:) ];
                ninds = zeros( size(newNodes,1), 1);
                ninds(i2)=i1;
                ninds(nidx)=noi;
                el2 = newElems;
                el2(:) = ninds( newElems(:) );
                obj.elems = [ obj.elems; el2 ];
          end
        end
        function el2 = mergeMesh( obj, newMesh )
            el2 = obj.merge(newMesh.nodes, newMesh.elems );
        end
        function el2 = append( obj, newNodes, newElems )
            nnodes=size(obj.nodes,1);
            obj.nodes = [ obj.nodes; newNodes ];
            obj.elems = [ obj.elems; newElems+nnodes ];
            
        end
        function el2 = connect( obj, sel, newNodes, newElems )
            [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
            newNodes = newNodes( si1, : );
            newElems = si2(newElems);
           if size( obj.nodes, 1 ) == 0 
                  obj.nodes = newNodes;
                  obj.elems = newElems;
                  el2 = newElems;
           else
                [~,i1,i2] = intersect( round(obj.nodes .* obj.tolerance), round(newNodes.*obj.tolerance), 'rows' );
                sn1 = find(sel.select( obj.nodes ));
                sn2 = find(sel.select( newNodes ));
                [C,ia] = setdiff( i2, sn2 );
                i1(ia) = [];
                i2(ia) = [];
                nidx  = 1:size(newNodes,1);
                nidx( i2 ) = [];
                noi   = 1:size(nidx,2);
                noi = noi + size(obj.nodes,1);
                obj.nodes = [ obj.nodes; newNodes(nidx,:) ];
                ninds = zeros( size(newNodes,1), 1);
                ninds(i2)=i1;
                ninds(nidx)=noi;
                el2 = newElems;
                el2(:) = ninds( newElems(:) );
                obj.elems = [ obj.elems; el2 ];
           end
            
        end
        function obj = transformDeg2D( obj, x0, angleDeg, xm )
            angleRad = angleDeg*pi/180;
            obj.nodes = [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                          x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad)] + xm;
                  
        end
        function nn = findClosestNode( obj, pt )
            nn = dsearchn(obj.nodes,pt);
        end
        function fnodes = findNodes(obj, selector)
            fnodes=find(selector.select(obj.nodes));
        end
        function felems = findElems( obj, selector, allElemNodes )
              fnodes = obj.findNodes( selector)';
              found = ismember(obj.elems, fnodes);
              if allElemNodes
                felems = find(sum(found,2)==size(obj.elems,2));
              else
                felems = find(sum(found,2)>0);
              end
        end

        function retelem = addRectMesh2D( obj, x1, y1, dx, dy, nx, ny, pattern )
            dim = max(max( pattern ));
            nn = ( dim *  nx + 1 ) * ( dim *  ny + 1 );
            ne = size( pattern, 1 ); 
            nfe = nx * ny;
            newNodes = zeros( nn, 2 );
            newElems = zeros( nfe, ne );
            newNodes(1:nn,1) = x1 + dx * rem( (1:nn) - 1, dim *  nx + 1  ) / dim / nx;
            newNodes(1:nn,2) = y1 + dy * ( floor( ( (1:nn) - 1 ) / ( dim *  nx + 1 ) ) / dim/ ny );
            gr = 1:nn;
            grid = reshape( gr, dim * nx + 1, dim * ny + 1 );
            for k=1:nfe
              ix = dim * rem((k)-1,nx)+1;
              iy = dim * floor(((k)-1)/(nx))+1;
              for i=1:ne
                newElems(k,i) = grid( ix + pattern(i,1),  iy + pattern(i,2) ) ;
              end
            end
            retelem = obj.merge( newNodes, newElems );
        end

        function obj = addDelaunayMesh2D( obj, P, C, nnodes )
            xv = unifrnd(min(P(:,1)),max(P(:,1)),1,nnodes);
            yv = unifrnd(min(P(:,2)),max(P(:,2)),1,nnodes);
            in = inpolygon(xv,yv,P(:,1),P(:,2));
            new_nodes=[xv(in); yv(in)]';
            new_nodes=[new_nodes; P];
            DT = delaunayTriangulation(new_nodes);
            ic = incenter(DT);
            in = inpolygon(ic(:,1),ic(:,2),P(:,1),P(:,2));
            %new_elems=delaunay(DT.ConnectivityList);
            %IO = isInterior(DT);
            obj.append( DT.Points, DT.ConnectivityList(in,:) );
            %obj.append( new_nodes, new_elems );
        end
        function obj = addRectMeshTriangular2D( obj, mode, x1, y1, dx, dy, nx, ny )
            ddx = dx/nx;
            ddy = dy/ny;
            offX =[0:ddx:(dx-ddx)] + x1;
            offY =[0:ddy:(dy-ddy)] + y1;
            [X,Y] = meshgrid(offX,offY);
            X=X';
            Y=Y';
            switch mode
                case 'quad'
                    newNodes = reshape([ 0+X(:) 0+Y(:) ddx+X(:) 0+Y(:) ddx/2+X(:) ddy/2+Y(:) ddx+X(:) ddy+Y(:) 0+X(:) ddy+Y(:) ]',2, 5*nx*ny)';
                    offN=[0:5:size(newNodes,1)-5]';
                    newElems = reshape([ 1+offN 2+offN 3+offN 2+offN 4+offN 3+offN 4+offN 5+offN 3+offN 5+offN 1+offN 3+offN ]',3,4*nx*ny)';
                case 'dual'
                    newNodes = reshape([ 0+X(:) 0+Y(:) ddx+X(:) 0+Y(:) 0+X(:) ddy+Y(:) ddx+X(:) ddy+Y(:) ]',2, 4*nx*ny)';
                    offN=[0:4:size(newNodes,1)-4]';
                    newElems = reshape([ 1+offN 2+offN 4+offN 1+offN 4+offN 3+offN ]',3,2*nx*ny)';
                otherwise
                error('no triangular generator in mode :'+mode)    ;
            end
            obj.merge( newNodes, newElems );
        end
        function obj = addRectMeshTetrahedral3D( obj, mode, x, dx, nx )
            ddx = dx./nx;

            offX =[0:ddx(1):(dx(1)-ddx(1))] + x(1);
            offY =[0:ddx(2):(dx(2)-ddx(2))] + x(2);
            offZ =[0:ddx(3):(dx(3)-ddx(3))] + x(3);
            
            [X1,Y1,Z1] = meshgrid(offY,offX,offZ);
            X=Y1(:);
            Y=X1(:);
            Z=Z1(:);
            switch mode
                case '6T'
                    nds=[ 0+X(:) 0+Y(:)      0+Z(:) ddx(1)+X(:) 0+Y(:)      0+Z(:) 0+X(:) ddx(2)+Y(:)      0+Z(:) ddx(1)+X(:) ddx(1)+Y(:)      0+Z(:)...
                          0+X(:) 0+Y(:) ddx(3)+Z(:) ddx(1)+X(:) 0+Y(:) ddx(3)+Z(:) 0+X(:) ddx(2)+Y(:) ddx(3)+Z(:) ddx(1)+X(:) ddx(1)+Y(:) ddx(3)+Z(:) ];
                    newNodes = reshape(nds',3, 8*nx(1)*nx(2)*nx(3))';
                    offN=[0:8:size(newNodes,1)-8]';
                    newElems = reshape([ 3+offN 7+offN 6+offN 5+offN ...
                                         3+offN 1+offN 5+offN 6+offN ...
                                         2+offN 1+offN 3+offN 6+offN ...
                                         6+offN 3+offN 7+offN 8+offN ...
                                         2+offN 6+offN 8+offN 3+offN ...
                                         2+offN 3+offN 8+offN 4+offN ]',4,6*nx(1)*nx(2)*nx(3))';
                otherwise
                error('no triangular generator in mode :'+mode)    ;
            end
            obj.merge( newNodes, newElems );
        end
        function newElems = addRectMesh3D( obj, x1, y1, z1, dx, dy, dz, nx, ny, nz, lnodes )
            ddx=dx/nx;
            ddy=dy/ny;
            ddz=dz/nz;
            nelems=nx*ny*nz;
            nnodes=size(lnodes,1)*nelems;
            newNodes=zeros(nnodes,size(lnodes,2));
            newElems=zeros(size(lnodes,1),nelems);
            newElems(:)=(1:nnodes);
            newElems = newElems';
            counter=1;
            for iz=1:nz
                for iy=1:ny
                    for ix=1:nx
                        newNodes(newElems(counter,:),1) = (2*x1+ddx*lnodes(:,1)+2*ddx*ix-ddx)/2;
                        newNodes(newElems(counter,:),2) = (2*y1+ddy*lnodes(:,2)+2*ddy*iy-ddy)/2;
                        newNodes(newElems(counter,:),3) = (2*z1+ddz*lnodes(:,3)+2*ddz*iz-ddz)/2;
                        counter=counter+1;
                    end
                end
            end
            obj.merge(newNodes, newElems);
        end

        function [newElems, zElems] = addRectMeshZlayer3D( obj, x1, y1, z1, dx, dy, dz, nx, ny, nz, zLayers, lnodes )
            ddx=dx/nx;
            ddy=dy/ny;
            ddz=dz/nz;
            nelems=nx*ny*nz;
            nnodes=size(lnodes,1)*nelems;
            newNodes=zeros(nnodes,size(lnodes,2));
            newElems=zeros(size(lnodes,1),nelems);
            newElems(:)=(1:nnodes);
            newElems = newElems';
            counter=1;
            zCounter=1;
            zElems=zeros(nx*ny,1);
            for iz=1:nz
                for iy=1:ny
                    for ix=1:nx
                        newNodes(newElems(counter,:),1) = (2*x1+ddx*lnodes(:,1)+2*ddx*ix-ddx)/2;
                        newNodes(newElems(counter,:),2) = (2*y1+ddy*lnodes(:,2)+2*ddy*iy-ddy)/2;
                        newNodes(newElems(counter,:),3) = (2*z1+ddz*lnodes(:,3)+2*ddz*iz-ddz)/2;
                        counter=counter+1;
                    end
                end
                if ( nz-iz < zLayers)
                   zElems(zCounter) = counter; 
                   zCounter=zCounter+1;
                end
            end
            zElems = obj.merge(newNodes, newElems);
        end

        function obj = addShapedMesh2D( obj, sf, x, ndiv1, ndiv2, pattern )
            baseMesh = Mesh();
            baseMesh.addRectMesh2D( -1, -1, 2, 2, ndiv1, ndiv2, pattern );
            baseMesh.nodes = sf.computeValue( baseMesh.nodes ) * x;
            obj.merge( baseMesh.nodes, baseMesh.elems );
        end
        function obj = addObjectMesh2D( obj, so, ndiv1, ndiv2, pattern )
            baseMesh = Mesh();
            baseMesh.addRectMesh2D( -1, -1, 2, 2, ndiv1, ndiv2, pattern );
            baseMesh.nodes = so.computeValue( baseMesh.nodes );
            obj.merge( baseMesh.nodes, baseMesh.elems );
        end
        function obj = addObjectMesh3D( obj, so, ndiv1, ndiv2, ndiv3, lnodes )
            baseMesh = Mesh();
            baseMesh.addRectMesh3D( -1, -1, -1, 2, 2, 2, ndiv1, ndiv2, ndiv3, lnodes );
            baseMesh.nodes = so.computeValue( baseMesh.nodes );
            obj.merge( baseMesh.nodes, baseMesh.elems );
        end
        function obj = addShapedMesh3D( obj, sf, x, ndiv, pattern )
            baseMesh = Mesh();
            baseMesh.addRectMesh3D( -1, -1, -1, 2, 2, 2, ndiv(1), ndiv(2), ndiv(3), pattern );
            baseMesh.nodes = sf.computeValue( baseMesh.nodes ) * x;
            obj.merge( baseMesh.nodes, baseMesh.elems );
        end
        function obj = addRectWithHoleInCornerMesh2D( obj, ri, x00, y00, s, div, pattern )
            xs = x00 - s;
            ys = y00 + s;
            alpha = 240*pi/180; 
            beta  = 315*pi/180;
            x1 = xs+ri*cos(alpha);
            mesh = Mesh();
            mesh.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            mesh.nodes = [ 0.5*( 1-mesh.nodes(:,2) ).*((x00-x1)/2*mesh.nodes(:,1)+(x00+x1)/2)+0.5*(mesh.nodes(:,2)+1).*(xs+ri*cos( (beta-alpha)/2 .* mesh.nodes(:,1) + (beta+alpha)/2)) ...
                           0.5*( 1-mesh.nodes(:,2) ).*(y00)+0.5*(mesh.nodes(:,2)+1).*(ys+ri*sin( (beta-alpha)/2 .* mesh.nodes(:,1) + (beta+alpha)/2)) ];
            obj.merge( mesh.nodes, mesh.elems );
            mesh = Mesh();
            mesh.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            alpha1 = alpha;           
            alpha = -45*pi/180;
            beta  = pi/6;         
            y1 = ys+ri*sin(beta);
            mesh.nodes = [ 0.5*( mesh.nodes(:,1)+1 ).*x00+0.5*(1-mesh.nodes(:,1)).*(xs+ri*cos( (beta-alpha)/2 .* mesh.nodes(:,2) + (beta+alpha)/2)) ...
                           0.5*( mesh.nodes(:,1)+1 ).*((y1-y00)/2*mesh.nodes(:,2)+(y00+y1)/2)+0.5*(1-mesh.nodes(:,1)).*(ys+ri*sin( (beta-alpha)/2 .* mesh.nodes(:,2) + (beta+alpha)/2)) ];           
            obj.merge( mesh.nodes, mesh.elems );
            obj.addShapedMesh2D( ShapeFunctionQ4(), [80 y00; x1 y00; 80 40; x1 ys+ri*sin(alpha1)], [div div], pattern );
            obj.addShapedMesh2D( ShapeFunctionQ4(),  [xs+ri*cos(beta) ys+ri*sin(beta); x00 y1; x00-40 x00-80; x00 x00-80], [div div], pattern );
        end
        function obj = addRectWithHoleMesh2D( obj, a, x0, y0, holefactor, div, pattern )
            r=a*holefactor;
            mesh = Mesh();
            mesh.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            mesh.nodes = [ 0.5*( 1-mesh.nodes(:,1)).*(x0-a)+0.5*(mesh.nodes(:,1)+1).*(x0-r*cos( pi/8 .* mesh.nodes(:,2) + pi/8)) ...
                           0.5*( 1-mesh.nodes(:,1)).*(a/2.*mesh.nodes(:,2)+y0+a/2)+0.5*(mesh.nodes(:,1)+1).*(y0+r.*sin( pi/8.*mesh.nodes(:,2)+pi/8)) ];
            mesh2 = Mesh();
            mesh2.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            mesh2.nodes  = [ 0.5*( mesh2.nodes(:,1)+1).*(x0+a)+0.5*(1-mesh2.nodes(:,1)).*(x0+r*cos( pi/8 .* mesh2.nodes(:,2) + pi/8)) ...
                         0.5*( mesh2.nodes(:,1)+1).*(a/2.*mesh2.nodes(:,2)+y0+a/2)+0.5*(1-mesh2.nodes(:,1)).*(y0+r.*sin( pi/8.*mesh2.nodes(:,2)+pi/8)) ];

            mesh.merge( mesh2.nodes, mesh2.elems );                              
            mesh.duplicateTransformedMeshDeg2D( [x0 y0], 90, [0 0] );
            mesh.duplicateTransformedMeshDeg2D( [x0 y0], 180, [0 0] );
            
            obj.merge(mesh.nodes,mesh.elems);
        end
        function obj = addRing2D( obj, x0, y0 , r1, r2, nr, nfi, pattern )
             mesh = Mesh();
             mesh.addRectMesh2D( r1, 0, r2-r1, 2*pi, nr, nfi, pattern );
             newNodes = [ x0+mesh.nodes(:,1).*cos( mesh.nodes(:,2) ) y0+mesh.nodes(:,1).*sin( mesh.nodes(:,2) ) ];
             obj.merge( newNodes, mesh.elems );
        end
        function obj = addQuarterCircle( obj, x0 , R, nr, pattern )
             mesh1 = Mesh();
             shapeFn = ShapeFunctionQ4();
             mesh1.addShapedMesh2D( shapeFn, x0+[ 0 0; R/2 0; 0 R/2; R/2/1.41 R/2/1.41 ], nr, pattern );

             sfL2 = ShapeFunctionL2();
             sl1 = ShapeObjectRectangular(sfL2,[0 R/2; R/2/1.41 R/2/1.41 ]);
             sl2 = ShapeObjectRectangular(sfL2,[R/2/1.41 R/2/1.41; R/2 0 ]);
             sc1 = CircleObject(x0,R,90,45);
             sc2 = CircleObject(x0,R,45,0);

             ms1 = MorphSpace(sl1,sc1);
             ms2 = MorphSpace(sl2,sc2);
                
             mesh1.addObjectMesh2D( ms1, nr, nr, pattern );
             mesh1.addObjectMesh2D( ms2, nr, nr, pattern );
             %mesh1.addShapedMesh2D( shapeFn, x0+[ R/2 0; R 0; R/2/1.41 R/2/1.41; R/1.41 R/1.41 ], nr, pattern );
             %mesh1.addShapedMesh2D( shapeFn, x0+[ 0 R/2; R/2/1.41 R/2/1.41; 0 R; R/1.41 R/1.41  ], nr, pattern );
             obj.merge( mesh1.nodes, mesh1.elems );
        end
        function obj = addHalfCircle( obj, x0 , R, nr, pattern )
            obj.addQuarterCircle( x0 , R, nr, pattern );
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            obj.addQuarterCircle( x0 , R, nr, pattern );
        end
        function obj = addCircle( obj, x0 , R, nr, pattern )
            obj.addQuarterCircle( x0 , R, nr, pattern );
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            obj.addQuarterCircle( x0 , R, nr, pattern );
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            obj.addQuarterCircle( x0 , R, nr, pattern );
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            obj.addQuarterCircle( x0 , R, nr, pattern );
        end
        function obj = addQuarterCylinder(obj, x0, R, h, nr, pattern)
                mesh1 = Mesh();
            
                a = R/(2*sqrt(2));
            
                shapeFn = ShapeFunctionH8();
                mesh1.addShapedMesh3D( ...
                    shapeFn, ...
                    x0 + [ 0   0   0;
                           R/2 0   0;
                           0   R/2 0;
                           a   a   0;
                           0   0   h;
                           R/2 0   h;
                           0   R/2 h;
                           a   a   h], ...
                    [nr(1) nr(1) nr(2)], pattern );
            
                sfL2 = ShapeFunctionQ4();
                sl1 = ShapeObjectRectangular(sfL2, x0 + [0   R/2 0;
                                                        a   a   0;
                                                        0   R/2 h;
                                                        a   a   h]);
                sl2 = ShapeObjectRectangular(sfL2, x0 + [a   a   0;
                                                        R/2 0   0;
                                                        a   a   h;
                                                        R/2 0   h]);
            
                % --- IMPORTANT: use absolute z-range for cylinder ---
                z0 = x0(3);
                z1 = x0(3) + h;
            
                sc1 = CylinderObject(x0, R, 90, 45, z0, z1);
                sc2 = CylinderObject(x0, R, 45, 0,  z0, z1);
            
                ms1 = MorphSpace(sl1, sc1);
                ms2 = MorphSpace(sl2, sc2);
            
                mesh1.addObjectMesh3D(ms1, nr(1), nr(1), nr(2), pattern);
                mesh1.addObjectMesh3D(ms2, nr(1), nr(1), nr(2), pattern);
            
                obj.merge(mesh1.nodes, mesh1.elems);
            end

        function obj = addQuarterCylinderOld( obj, x0 , R, h, nr, pattern )
             mesh1 = Mesh();
             shapeFn = ShapeFunctionH8();
             mesh1.addShapedMesh3D( shapeFn, x0+[ 0 0 0; R/2 0 0; 0 R/2 0; R/2/1.41 R/2/1.41 0; 0 0 h; R/2 0 h; 0 R/2 h; R/2/1.41 R/2/1.41 h], [nr(1) nr(1) nr(2)], pattern );

             sfL2 = ShapeFunctionQ4();
             a = R/(2*sqrt(2));

            sl1 = ShapeObjectRectangular(sfL2, x0+[0 R/2 0; a a 0; 0 R/2 h; a a h]);
            sl2 = ShapeObjectRectangular(sfL2, x0+[a a 0; R/2 0 0; a a h; R/2 0 h]);

           
             sc1 = CylinderObject(x0,R,90,45,0,h);
             sc2 = CylinderObject(x0,R,45,0,0,h);

             ms1 = MorphSpace(sl1,sc1);
             ms2 = MorphSpace(sl2,sc2);
                
             mesh1.addObjectMesh3D( ms1, nr(1), nr(1), nr(2), pattern );
             mesh1.addObjectMesh3D( ms2, nr(1), nr(1), nr(2), pattern );
             obj.merge( mesh1.nodes, mesh1.elems );
        end

        function obj = addLayeredQuarterCylinder(obj, x0, R, h, nrZ, nrXY, add_interface, i_th, pattern)
            %ADDLAYEREDQUARTERCYLINDER Stack quarter-cylinder layers in Z with optional interface layers.
            %
            %   obj = addLayeredQuarterCylinder(obj, x0, R, h, nrZ, nrXY, add_interface, i_th, pattern)
            %
            % Inputs:
            %   x0            - base origin [x y z] of the cylinder stack
            %   R             - radius
            %   h             - vector of nominal layer heights (sum(h) = total height)
            %   nrZ           - vector of z-resolutions per nominal layer (same length as h)
            %   nrXY          - scalar in-plane resolution (used as nr(1) in addQuarterCylinder)
            %   add_interface - logical; if true, insert interface layers between nominal layers
            %   i_th          - interface thickness (height), used only if add_interface==true
            %   pattern       - pattern passed to mesh generation
            %
            % Behavior with add_interface==true:
            %   - Inserts (L-1) interface layers of height i_th, each with nrZ_interface = 1.
            %   - To keep total height unchanged (still sum(h)), each nominal layer k is shortened:
            %       h_eff(k) = h(k) - (k>1)*i_th/2 - (k<L)*i_th/2
            %     i.e., remove i_th/2 at each internal boundary from the adjacent layers.
            %
            % Notes:
            %   - z coordinates are stacked from x0(3) upwards.
            %   - Uses existing obj.addQuarterCylinder(x0, R, h_layer, [nrXY nrZ_layer], pattern)
            
                % --- sanity checks ---
                if numel(x0) ~= 3
                    error('x0 must be a 1x3 vector.');
                end
                if ~isvector(h) || ~isvector(nrZ)
                    error('h and nrZ must be vectors.');
                end
            
                h   = h(:);
                nrZ = nrZ(:);
            
                L = numel(h);
                if numel(nrZ) ~= L
                    error('h and nrZ must have the same number of layers.');
                end
                if any(h <= 0)
                    error('All layer heights h must be > 0.');
                end
                if any(nrZ < 1) || any(mod(nrZ,1) ~= 0)
                    error('All nrZ entries must be positive integers.');
                end
                if ~(isscalar(nrXY) && nrXY >= 1 && mod(nrXY,1) == 0)
                    error('nrXY must be a positive integer scalar.');
                end
                if ~(islogical(add_interface) || (isscalar(add_interface) && (add_interface==0 || add_interface==1)))
                    error('add_interface must be logical true/false.');
                end
                add_interface = logical(add_interface);
            
                if add_interface
                    if ~(isscalar(i_th) && i_th > 0)
                        error('i_th must be a positive scalar when add_interface is true.');
                    end
                    if L < 2
                        % No interfaces possible; just build the single layer normally
                        add_interface = false;
                    end
                end
            
                % --- compute effective layer heights if interfaces are enabled ---
                if add_interface
                    h_eff = h;
                    for k = 1:L
                        if k > 1
                            h_eff(k) = h_eff(k) - i_th/2;
                        end
                        if k < L
                            h_eff(k) = h_eff(k) - i_th/2;
                        end
                    end
            
                    if any(h_eff <= 0)
                        error('Interface thickness i_th is too large: some effective layer heights become <= 0.');
                    end
                else
                    h_eff = h;
                end
            
                % --- build stack ---
                zoff = 0.0;
            
                for k = 1:L
                    % main layer k
                    obj = obj.addQuarterCylinder( x0 + [0 0 zoff], R, h_eff(k), [nrXY nrZ(k)], pattern );
                    zoff = zoff + h_eff(k);
            
                    % interface between k and k+1
                    if add_interface && (k < L)
                        obj = obj.addQuarterCylinder( x0 + [0 0 zoff], R, i_th, [nrXY 1], pattern );
                        zoff = zoff + i_th;
                    end
                end
            
                % --- optional consistency check (floating tolerance) ---
                targetH = sum(h);
                if abs(zoff - targetH) > 1e-10 * max(1, targetH)
                    error('Height mismatch: built height %.15g, expected %.15g.', zoff, targetH);
                end
            end


        function obj = addLayeredQuarterCylinderInt(obj, x0, R, h, nrZ, nrXY, add_interface, i_th, pattern)
                zoff = 0;
                for k = 1:numel(h)
                    obj = obj.addQuarterCylinder( x0 + [0 0 zoff], R, h(k), [nrXY nrZ(k)], pattern );
                    zoff = zoff + h(k);
                end
                %obj.merge(mesh1.nodes, mesh1.elems);
        end



        function obj = addHalfCylinder( obj, x0 , R, h, nr, localNodes )
            obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
        end
        function obj = addCylinder( obj, x0 , R, h, nr, localNodes )
            obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
        end
        function addPipe3D(obj,x0,r,R,al1,al2,h1,h2,nr,nc,nz,lnodes)
            mesh=Mesh();
            mesh.addRectMesh3D( r, deg2rad(al1), h1, R-r, deg2rad(al2-al1), h2-h1, nr, nc, nz, lnodes);
            mesh.transformToCylindrical3D(x0);
            obj.mergeMesh(mesh);
        end

        function addManipulatorHalfSegment3D(obj,r,R,h,alpha,nr,nc,nz,lnodes)
            mesh=Mesh();
            mesh.addRectMesh3D( r, 0, 0, R-r, 2*pi, h, nr, nc, nz, lnodes);
            mesh.transformToCylindrical3D([R 0]);
           [mesh.convex_hull_nodes, ~] = convhull(mesh.nodes(:,1), mesh.nodes(:,2), mesh.nodes(:,3));

            xp=mesh.nodes;
            xp(:,3)=h;
            c=cos(alpha);
            s=sin(alpha);
            t=tan(alpha);
            Ry=[c 0 s; 0 1 0; -s 0 c];
            nnodes=[0 0 h]+(Ry*(xp-[0 0 h])')';
            %x=(mesh.nodes(:,3)/h).*nnodes(:,1)+(1-mesh.nodes(:,3)/h).*mesh.nodes(:,1);
            x=nnodes(:,1);
            z=(mesh.nodes(:,3)/h).*(h-mesh.nodes(:,1)*t);
            mesh.nodes=[x mesh.nodes(:,2) z];
            obj.mergeMesh(mesh);
            obj.convex_hull_nodes=mesh.convex_hull_nodes;
        end

        function obj = addrectPipe(obj,w,h,l1,th,nth,lnodes)
            nx=round(l1/th)*nth;
            ny=round((w-2*th)/th)*nth;
            nz=round((h-2*th)/th)*nth;
            obj.addRectMesh3D( 0, 0, 0, l1, th, th, nx, nth, nth, lnodes );
            obj.addRectMesh3D( 0, 0, h-th, l1, th, th, nx, nth, nth, lnodes );
            obj.addRectMesh3D( 0, w-th, h-th, l1, th, th, nx, nth, nth, lnodes );
            obj.addRectMesh3D( 0, w-th, 0, l1, th, th, nx, nth, nth, lnodes );

            obj.addRectMesh3D( 0, 0, th, l1, th, h-2*th, nx, nth, nz, lnodes );
            obj.addRectMesh3D( 0, w-th, th, l1, th, h-2*th, nx, nth, nz, lnodes );

            obj.addRectMesh3D( 0, th, 0, l1, w-2*th, th, nx, ny, nth, lnodes );
            obj.addRectMesh3D( 0, th, h-th, l1, w-2*th, th, nx, ny, nth, lnodes );
        end

        function obj = addLshape( obj, l, h, nh, pattern )
            nl = round(l/h*nh+0.5);
            obj.addRectMesh2D( 0, 0, h, h, nh, nh, pattern );
            obj.addRectMesh2D( 0, h, h, l-h, nh, nl-nh, pattern );
            obj.addRectMesh2D( h, 0, l-h, h, nl-nh, nh, pattern );
        end

        function obj = addLshape3D( obj, l, h, nh, lnodes )
            nl = round(l/h*nh+0.5);
            obj.addRectMesh3D( 0, 0, 0, h,   h, h,   nh,    nh, nh,    lnodes);
            obj.addRectMesh3D( 0, 0, h, h,   h, l-h, nh,    nh, nl-nh, lnodes);
            obj.addRectMesh3D( h, 0, 0, l-h, h, h,   nl-nh, nh, nh,    lnodes );
        end
        function obj = duplicateTransformedMeshDeg2D( obj, x0, angleDeg, xm )
            angleRad = angleDeg*pi/180;
            obj.merge( [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                         x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad)] + xm, obj.elems );
                  
        end
        function obj = duplicateTransformedMeshDeg3D( obj, x0, angleDeg, xm )
            angleRad = angleDeg*pi/180;
            obj.merge( [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                         x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad)...
                         obj.nodes(:,3)] + xm, obj.elems );
                  
        end
        function obj = array( obj, coord, n )
                m = max(obj.nodes)-min(obj.nodes);
                mv=zeros(1,size(obj.nodes,2));
                belems=obj.elems;
                bnodes=obj.nodes;
                for k=1:n
                    mv(coord)=k*m(coord);
                    newnodes=bnodes+mv;
                    obj.merge( newnodes, belems );
                end

        end
        function obj = transformMeshDeg2D( obj, x0, angleDeg, xm )
            angleRad = angleDeg*pi/180;
            obj.nodes = [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                          x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad)] + xm;
                  
        end

        function cone_nodes = coneTransformationX(obj, l, r1, r2, nodes)
            %CONETRANSFORMATIONX Taper (cone) transform around Z axis.
            %
            %   cone_nodes = coneTransformationX(l, r1, r2, nodes)
            %
            % Inputs:
            %   l     - height of cone (z from 0 to l)
            %   r1    - radius at z = 0
            %   r2    - radius at z = l
            %   nodes - Nx3 array of [x y z]
            %
            % Output:
            %   cone_nodes - Nx3 transformed nodes
            %
            % Behavior:
            %   - z is unchanged
            %   - x,y are scaled by s(z) = 1 + (r2/r1 - 1) * (z/l)
            %     so at z=0:  s=1
            %        z=l:  s=r2/r1  (radius changes from r1 -> r2)
            
                cone_nodes = nodes;
            
                % Basic checks
                if size(nodes,2) ~= 3
                    error('nodes must be an Nx3 array.');
                end
                if l <= 0
                    error('l must be > 0.');
                end
                if r1 == 0
                    error('r1 must be non-zero (cannot scale by r2/r1).');
                end
            
                z = nodes(:,3);
            
                % Clamp interpolation to [0, 1] in case some nodes are slightly outside
                t = z ./ l;
                t = max(0, min(1, t));
            
                % Linear scale from 1 at z=0 to (r2/r1) at z=l
                s = 1 + (r2/r1 - 1) .* t;
            
                % Apply scaling in XY plane
                cone_nodes(:,1) = nodes(:,1) .* s;
                cone_nodes(:,2) = nodes(:,2) .* s;
        end

        function obj = transformMesh3DDegXY( obj, x0, angleDeg, xm )
            angleRad = angleDeg*pi/180;
            obj.nodes = [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                          x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad) obj.nodes(:,3)] + xm;
                  
        end
        function obj = transformToPolar2D( obj, x0, y0 )
             newNodes = [ x0+obj.nodes(:,1).*cos( obj.nodes(:,2) ) y0+obj.nodes(:,1).*sin( obj.nodes(:,2) ) ];
             [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
            obj.nodes = newNodes( si1, : );
            obj.elems = si2(obj.elems);
        end
        function obj = transformToCylindrical3D( obj, x0 )
             newNodes = [ x0(1)+obj.nodes(:,1).*cos( obj.nodes(:,2) ) x0(2)+obj.nodes(:,1).*sin( obj.nodes(:,2) ) obj.nodes(:,3) ];
             [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
            obj.nodes = newNodes( si1, : );
            obj.elems = si2(obj.elems);
        end
        function removeNodes( obj, selector )
            nodesToRemove=find(selector.select(obj.nodes));
            nodesToLeave=find(1-selector.select(obj.nodes));
            elemsToRemove=find(sum(ismember(obj.elems,nodesToRemove),2)>0);
            oldnn=size(obj.nodes,1);
            numbers=zeros(oldnn,1);
            numbers(nodesToLeave)=1:size(nodesToLeave,1);
            obj.elems=numbers(obj.elems);
            obj.nodes(nodesToRemove,:)=[];
            obj.elems(elemsToRemove,:)=[];
        end
        function removeNodesByNumbers( obj, nodesToRemove )
            nodesToLeave=1:size(obj.nodes,1);
            nodesToLeave(nodesToRemove)=[];
            elemsToRemove=find(sum(ismember(obj.elems,nodesToRemove),2)>0);
            oldnn=size(obj.nodes,1);
            numbers=zeros(oldnn,1);
            numbers(nodesToLeave)=1:size(nodesToLeave,1);
            obj.elems=numbers(obj.elems);
            obj.nodes(nodesToRemove,:)=[];
            obj.elems(elemsToRemove,:)=[];
        end
        function removeElemsByNumbers( obj, elemsToRemove )
            elemsToLeave=(1:size(obj.elems,1))';
            elemsToLeave(elemsToRemove)=[];
            nrem=obj.elems(elemsToLeave,:);
            nodesToLeave=unique(nrem(:));
            nodesToRemove=(1:size(obj.nodes,1))';
            nodesToRemove(nodesToLeave)=[];
            oldnn=size(obj.nodes,1);
            mapNodeNumbers=zeros(oldnn,1);
            mapNodeNumbers(nodesToLeave)=1:size(nodesToLeave,1);
            obj.elems=mapNodeNumbers(obj.elems);
            obj.nodes(nodesToRemove,:)=[];
            obj.elems(elemsToRemove,:)=[];
        end
        function transformNodesXY( obj, transformFn )
            obj.nodes = transformFn( obj.nodes );
        end
        function laplacianSmoothing(n)
            
        end
        function importFEMesh( obj, mesh )
            obj.nodes=mesh.Nodes';
            obj.elems=mesh.Elements';
        end
        function exportMeshToFile(obj, filenamebase)
            nodes=obj.nodes;
            elems=obj.elems;
            save([filenamebase '_mesh.mat'],"nodes","elems");
            dlmwrite([filenamebase '_nodes.txt'],obj.nodes,'delimiter','\t','precision','%7.3f');
            dlmwrite([filenamebase '_elems.txt'],obj.elems,'delimiter','\t','precision',6,'-append');
        end
        function upward_facing_nodes = findUpwardFacingNodes(obj)
            % Preallocate for storing upward-facing node indices
            upward_facing_nodes = [];
            
            % Define a threshold for the Z-component of the normal
            z_threshold = 0.1;  % Adjust this value as needed (e.g., 0.1 means the normal must be at least 10% upward)

            % Compute the convex hull
            k=obj.convex_hull_nodes;
    
            % Loop over each facet to calculate the normals and check if they are upward-facing
            for i = 1:size(k,1)
                % Get the indices of the vertices for the i-th facet
                idx = k(i,:);
                
                % Get the vertices of the triangle
                v1 = obj.nodes(idx(1),:);
                v2 = obj.nodes(idx(2),:);
                v3 = obj.nodes(idx(3),:);
                
                % Compute edge vectors
                edge1 = v2 - v1;
                edge2 = v3 - v1;
                
                % Compute the normal using the cross product
                normal = cross(edge1, edge2);
                
                % Normalize the normal vector
                normal = normal / norm(normal);
                
                % Check if the normal vector is pointing upwards (Z-component above threshold)
                if normal(3) > z_threshold
                    % Add the vertices of this facet to the upward-facing nodes list
                    upward_facing_nodes = [upward_facing_nodes; idx'];
                end
            end
            
            % Remove duplicate nodes (since a node can belong to multiple upward-facing facets)
            upward_facing_nodes = unique(upward_facing_nodes);
        end
        function downward_facing_nodes = findDownwardFacingNodes(obj)
            % Preallocate for storing upward-facing node indices
            downward_facing_nodes = [];
            
            % Define a threshold for the Z-component of the normal
            z_threshold = 0.1;  % Adjust this value as needed (e.g., 0.1 means the normal must be at least 10% upward)

            % Compute the convex hull
            k=obj.convex_hull_nodes;
    
            % Loop over each facet to calculate the normals and check if they are upward-facing
            for i = 1:size(k,1)
                % Get the indices of the vertices for the i-th facet
                idx = k(i,:);
                
                % Get the vertices of the triangle
                v1 = obj.nodes(idx(1),:);
                v2 = obj.nodes(idx(2),:);
                v3 = obj.nodes(idx(3),:);
                
                % Compute edge vectors
                edge1 = v2 - v1;
                edge2 = v3 - v1;
                
                % Compute the normal using the cross product
                normal = cross(edge1, edge2);
                
                % Normalize the normal vector
                normal = normal / norm(normal);
                
                % Check if the normal vector is pointing upwards (Z-component above threshold)
                if normal(3) < -z_threshold
                    % Add the vertices of this facet to the upward-facing nodes list
                    downward_facing_nodes = [downward_facing_nodes; idx'];
                end
            end
            
            % Remove duplicate nodes (since a node can belong to multiple upward-facing facets)
            downward_facing_nodes = unique(downward_facing_nodes);
        end
    end
end

