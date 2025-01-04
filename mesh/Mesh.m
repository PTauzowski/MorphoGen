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
            obj.addShapedMesh2D( ShapeFunctionL4(), [80 y00; x1 y00; 80 40; x1 ys+ri*sin(alpha1)], [div div], pattern );
            obj.addShapedMesh2D( ShapeFunctionL4(),  [xs+ri*cos(beta) ys+ri*sin(beta); x00 y1; x00-40 x00-80; x00 x00-80], [div div], pattern );
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
             shapeFn = ShapeFunctionL4();
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
        function obj = addQuarterCylinder( obj, x0 , R, h, nr, pattern )
             mesh1 = Mesh();
             shapeFn = ShapeFunctionL8();
             mesh1.addShapedMesh3D( shapeFn, x0+[ 0 0 0; R/2 0 0; 0 R/2 0; R/2/1.41 R/2/1.41 0; 0 0 h; R/2 0 h; 0 R/2 h; R/2/1.41 R/2/1.41 h], [nr(1) nr(1) nr(2)], pattern );

             sfL2 = ShapeFunctionL4();
             sl1 = ShapeObjectRectangular(sfL2,x0+[0 R/2 0; R/2/1.41 R/2/1.41 0; 0 R/2 h; R/2/1.41 R/2/1.41 h ]);
             sl2 = ShapeObjectRectangular(sfL2,x0+[R/2/1.41 R/2/1.41 0; R/2 0 0; R/2/1.41 R/2/1.41 h; R/2 0 h]);
             sc1 = CylinderObject(x0,R,90,45,0,h);
             sc2 = CylinderObject(x0,R,45,0,0,h);

             ms1 = MorphSpace(sl1,sc1);
             ms2 = MorphSpace(sl2,sc2);
                
             mesh1.addObjectMesh3D( ms1, nr(1), nr(1), nr(2), pattern );
             mesh1.addObjectMesh3D( ms2, nr(1), nr(1), nr(2), pattern );
             %mesh1.addShapedMesh2D( shapeFn, x0+[ R/2 0; R 0; R/2/1.41 R/2/1.41; R/1.41 R/1.41 ], nr, pattern );
             %mesh1.addShapedMesh2D( shapeFn, x0+[ 0 R/2; R/2/1.41 R/2/1.41; 0 R; R/1.41 R/1.41  ], nr, pattern );
             obj.merge( mesh1.nodes, mesh1.elems );
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
            mesh.addRectMesh3D( r, deg2rad(al1), h1, R-r, deg2rad(al2-al1), h2-h1, nr, nc, nz, lnodes)
            mesh.transformToCylindrical3D(x0);
            obj.mergeMesh(mesh);
        end

        function addManipulatorHalfSegment3D(obj,r,R,h,alpha,nr,nc,nz,lnodes)
            mesh=Mesh();
            mesh.addRectMesh3D( r, 0, 0, R-r, 2*pi, h, nr, nc, nz, lnodes);
            mesh.transformToCylindrical3D([0 0]);
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
        function leaveElemsByNumbers( obj, elemsToLeave )
            elemsToRemove=true(size(obj.elems,1),1);
            elemsToRemove(elemsToLeave)=false;
            obj.removeElemsByNumbers(elemsToRemove);
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
        
        % function [tnodes, telems] = getTetrahedralMesh(obj,selected)
        %     Tetrahedra = [
        %         2 5 3 1;  % Tetrahedron 1
        %         2 8 3 5;  % Tetrahedron 2
        %         3 5 8 7;  % Tetrahedron 3
        %         8 5 2 6;  % Tetrahedron 4
        %         3 8 2 4;  % Tetrahedron 5
        %     ];
        %     snodes=false(size(obj.nodes,1),1);
        %     snodes(obj.elems(selected,:))=true;
        %     nSelElems=size(find(selected),1);
        %     nSelNodes=size(find(snodes),1);
        %     tnodes=obj.nodes;
        %     telems=reshape(obj.elems(selected,Tetrahedra'),4,5*nSelElems)';
        % end

        function [tnodes, telems] = getTetrahedralMesh(obj, selected)
            % Tetrahedra mapping for a single hexahedron
            Tetrahedra = [
                5 2 1 4;  % Tetrahedron 1
                8 5 2 6;  % Tetrahedron 2
                4 5 2 8;  % Tetrahedron 3
                7 5 4 8;  % Tetrahedron 4
                4 7 1 5;  % Tetrahedron 5
                1 7 4 3;  % Tetrahedron 5
            ];
            
            % Identify selected nodes
            snodes = false(size(obj.nodes, 1), 1);
            snodes(obj.elems(selected, :)) = true;
            
            % Map global indices of selected nodes
            globalToLocal = find(snodes);
            localToGlobal = zeros(size(snodes));
            localToGlobal(globalToLocal) = 1:length(globalToLocal);
            
            % Filter nodes and elements
            tnodes = obj.nodes(globalToLocal, :);
            nSelElems = sum(selected);
            
            % Generate tetrahedrons for each selected hexahedron
            telems = zeros(nSelElems * size(Tetrahedra, 1), 4); % Each hexahedron creates 5 tetrahedrons
            cnt = 1;
            for i = find(selected)'
                hexNodes = obj.elems(i, :); % Global node indices of current hexahedron
                localHexNodes = localToGlobal(hexNodes); % Map to local indices
                telems(cnt:cnt + 5, :) = localHexNodes(Tetrahedra);
                cnt = cnt + 6;
            end
        end

        function exportMeshToFile(obj, selected, filenamebase)
            %save([filenamebase '_mesh.mat'],"nodes","elems");
            dlmwrite(filenamebase + '_hexamesh.txt',obj.nodes,'delimiter',',','precision','%7.3f');
            dlmwrite(filenamebase + '_hexamesh.txt',obj.elems(selected,:),'delimiter','\t','precision',6,'-append');
            
            [~, telems] = obj.getTetrahedralMesh(selected);
            
            dlmwrite(filenamebase + '_tetramesh.txt',obj.nodes,'delimiter',',','precision','%7.3f');
            dlmwrite(filenamebase + '_tetramesh.txt',telems,'delimiter',',','precision',6,'-append');
        end
        function exportToPLY(obj,filename)
            % Open the file
            fid = fopen(filename+".ply", 'w');
            if fid == -1
                error('Could not create file');
            end
            
            % Write PLY header
            fprintf(fid, 'ply\nformat ascii 1.0\n');
            fprintf(fid, 'element vertex %d\n', size(obj.nodes, 1));
            fprintf(fid, 'property float x\nproperty float y\nproperty float z\n');
            fprintf(fid, 'element face %d\n', size(obj.elems, 1) * 4); % 4 faces per tetrahedron
            fprintf(fid, 'property list uchar int vertex_indices\n');
            fprintf(fid, 'end_header\n');
            
            % Write vertices
            fprintf(fid, '%.6f %.6f %.6f\n', obj.nodes');
            
            % Write faces (export tetrahedrons as triangular faces)
            for i = 1:size(obj.elems, 1)
                fprintf(fid, '3 %d %d %d\n', obj.elems(i, [1, 2, 3]) - 1); % Face 1
                fprintf(fid, '3 %d %d %d\n', obj.elems(i, [1, 3, 4]) - 1); % Face 2
                fprintf(fid, '3 %d %d %d\n', obj.elems(i, [1, 4, 2]) - 1); % Face 3
                fprintf(fid, '3 %d %d %d\n', obj.elems(i, [2, 3, 4]) - 1); % Face 4
            end
            
            % Close the file
            fclose(fid);
            disp(['Mesh exported to ', filename]);
        end
        function exportTetraToSTL(obj,filename)
            % Open the file for writing
            fid = fopen(filename+".stl", 'w');
            if fid == -1
                error('Could not create file');
            end
            
            % Write STL header
            fprintf(fid, 'solid TetrahedralMesh\n');
            
            % Extract surface faces
            % Each tetrahedron face: 4 faces per tetrahedron
            tetraFaces = [
                1 2 3;  % Face 1
                1 3 4;  % Face 2
                1 4 2;  % Face 3
                2 3 4   % Face 4
            ];
            
            % Store unique surface faces
            faces = [];
            for i = 1:size(obj.elems, 1)
                tetra = obj.elems(i, :); % Extract tetrahedron nodes
                faces = [faces; tetra(tetraFaces)]; % Append all faces
            end
            
            % Remove duplicate faces (shared between neighboring tetrahedra)
            faces = unique(sort(faces, 2), 'rows');
            
            % Write each triangular face in STL format
            for i = 1:size(faces, 1)
                v1 = obj.nodes(faces(i, 1), :);
                v2 = obj.nodes(faces(i, 2), :);
                v3 = obj.nodes(faces(i, 3), :);
                
                % Compute normal vector for the triangle
                n = cross(v2 - v1, v3 - v1);
                n = n / norm(n);
                
                % Write the face
                fprintf(fid, '  facet normal %.6f %.6f %.6f\n', n);
                fprintf(fid, '    outer loop\n');
                fprintf(fid, '      vertex %.6f %.6f %.6f\n', v1);
                fprintf(fid, '      vertex %.6f %.6f %.6f\n', v2);
                fprintf(fid, '      vertex %.6f %.6f %.6f\n', v3);
                fprintf(fid, '    endloop\n');
                fprintf(fid, '  endfacet\n');
            end
            
            % Write STL footer
            fprintf(fid, 'endsolid TetrahedralMesh\n');
            
            % Close the file
            fclose(fid);
            disp(['Tetrahedral mesh surface exported to ', filename]);
        end
        function exportToStep(obj,filename)
            % EXPORTTOSTEP Exports a 3D mesh to a basic STEP file.
            % 
            % INPUTS:
            %   filename - Name of the STEP file to write (e.g., 'output.step').
            %   nodes - Matrix of node coordinates (Nx3).
            %   elems - Matrix of triangular elements (Mx3), referencing node indices.
            %
            % NOTE: This function generates a simplified STEP file with only B-Rep data.
            
            fid = fopen(filename+".step", 'w');
            if fid == -1
                error('Could not open file for writing.');
            end
        
            % STEP Header
            fprintf(fid, "ISO-10303-21;\n");
            fprintf(fid, 'HEADER;\n');
            fprintf(fid, "FILE_DESCRIPTION(('Basic Mesh Export'),'2;1');\n");
            fprintf(fid, "FILE_NAME('%s','%s');\n', filename, datestr(now, 'yyyy-mm-ddTHH:MM:SS\n");
            fprintf(fid, "FILE_SCHEMA(('AUTOMOTIVE_DESIGN_CC2'));\n");
            fprintf(fid, 'ENDSEC;\n');
            
            % STEP Data Section
            fprintf(fid, 'DATA;\n');
            
            % Write nodes as STEP vertices
            for i = 1:size(obj.nodes, 1)
                fprintf(fid, "#%d = CARTESIAN_POINT ( '' , ( %.6f, %.6f, %.6f ) );\n", ...
                    i, obj.nodes(i, 1), obj.nodes(i, 2), obj.nodes(i, 3));
            end
            
            % Write triangles as STEP faces
            faceId = size(obj.nodes, 1) + 1;
            for i = 1:size(obj.elems, 1)
                fprintf(fid, "#%d = ADVANCED_FACE ( '' , ( ", faceId + i - 1);
                for j = 1:3
                    edgeId = faceId + size(obj.elems, 1) + ((i - 1) * 3) + j;
                    fprintf(fid, '#%d ', edgeId);
                    if j < 3
                        fprintf(fid, ', ');
                    else
                        fprintf(fid, ')');
                    end
                end
                fprintf(fid, ', #4) ;\n');
            end
            
            % Close STEP file
            fprintf(fid, 'ENDSEC;\n');
            fprintf(fid, 'END-ISO-10303-21;\n');
            fclose(fid);
            
            disp(['STEP file exported to ', filename]);
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

