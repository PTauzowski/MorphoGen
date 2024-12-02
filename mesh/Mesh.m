classdef Mesh < handle
    
    properties
        nodes;
        tolerance = 10000;
    end
    
    methods
        function dim = getDim(obj)
            dim=size(obj.nodes,2);
        end
        function nn = getNumberOfNodes(obj)
            nn = size(obj.nodes,1);
        end
        function mergedElems = merge( obj, newNodes, newElems )
            [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
            newNodes = newNodes( si1, : );
            newElems = si2(newElems);
            if  size(newElems,2)==1
                newElems=newElems';
            end
           if size( obj.nodes, 1 ) == 0 
                  obj.nodes = newNodes;
                  mergedElems = newElems;
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
                mergedElems = newElems;
                mergedElems(:) = ninds( newElems(:) );
          end
        end
       
        function newIndices = mergeMesh( obj, newMesh )
            [~, iNodes, iNewNodes] = intersect( round(obj.nodes .* obj.tolerance), round(newMesh.nodes.*obj.tolerance), 'rows' );
            uniqueIndices=1:size(newMesh.nodes,1);
            uniqueIndices(iNewNodes)=[];
            newIndices=1:size(newMesh.nodes,1);
            newIndices(iNewNodes)=iNodes;
            newIndices(uniqueIndices)=1:size(uniqueIndices,2)+size(obj.nodes,1);
            obj.nodes=[obj.nodes; newMesh.nodes(uniqueIndices,:)];
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
        function felems = findElems( obj, elems, selector )
              fnodes = obj.findNodes( selector)';
              found = ismember(elems, fnodes);
              felems = find(sum(found,2)==size(elems,2));
        end
        function elems = addRectMesh2D( obj, x1, y1, dx, dy, nx, ny, pattern )
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
            elems = obj.merge( newNodes, newElems );
        end
        function elems = addRectMeshArray2D( obj, x1, y1, dx, dy, minres, pattern )
                elemsize = min([dx dy])/minres;
                xp=x1;
                yp=y1;
                mesh=Mesh();
                elems = [];
                for k=1:size(dx,2)
                    yp=y1;
                    for l=1:size(dy,2)
                        elems = [elems; mesh.addRectMesh2D(xp,yp,dx(k),dy(l),round(dx(k)/elemsize),round(dy(l)/elemsize),pattern)];
                        yp=yp+dy(l);
                    end
                    xp=xp+dx(k);
                end
                elems = obj.merge( mesh.nodes, elems );
        end
        function elems = addDelaunayMesh2D( obj, P, C, nnodes )
            xv = unifrnd(min(P(:,1)),max(P(:,1)),1,nnodes);
            yv = unifrnd(min(P(:,2)),max(P(:,2)),1,nnodes);
            in = inpolygon(xv,yv,P(:,1),P(:,2));
            new_nodes=[xv(in); yv(in)]';
            %new_nodes=[new_nodes; P];
            new_nodes=P;
            DT = delaunayTriangulation(new_nodes);
            ic = incenter(DT);
            in = inpolygon(ic(:,1),ic(:,2),P(:,1),P(:,2));
            %new_elems=delaunay(DT.ConnectivityList);
            IO = isInterior(DT);
            newNodes = obj.laplacianSmoothing(DT.Points, DT.ConnectivityList(in,:), 20);
            elems = obj.merge( newNodes, DT.ConnectivityList(in,:) );
            %obj.append( new_nodes, new_elems );
        end
        function smoothedNodes = laplacianSmoothing(obj, nodes, connectivityArray, numIterations)
            % nodes: N x 2 or N x 3 array of node coordinates
            % connectivityArray: M x 3 array of triangle vertex indices
            % numIterations: Number of smoothing iterations
            % smoothedNodes: the smoothed node coordinates
        
            % Number of nodes
            numNodes = size(nodes, 1);
            
            % Determine boundary nodes
            boundaryNodes = obj.detectBoundaryNodes(connectivityArray, numNodes);
            
            % Initialize smoothedNodes with the original nodes
            smoothedNodes = nodes;
        
            % Perform smoothing for the specified number of iterations
            for iter = 1:numIterations
                % Create arrays to accumulate neighbor sums for each dimension (x, y, z)
                neighborSumX = zeros(numNodes, 1);
                neighborSumY = zeros(numNodes, 1);
                
                % If 3D mesh, also create a z-component accumulator
                if size(nodes, 2) == 3
                    neighborSumZ = zeros(numNodes, 1);
                end
                
                neighborCount = zeros(numNodes, 1);
        
                % Get all edges from the connectivityArray
                edges = [connectivityArray(:, [1, 2]); connectivityArray(:, [2, 3]); connectivityArray(:, [3, 1])];
        
                % Accumulate neighbor sums and counts for X and Y (and Z if 3D)
                neighborSumX = accumarray([edges(:,1); edges(:,2)], [nodes(edges(:,2),1); nodes(edges(:,1),1)], [numNodes, 1]);
                neighborSumY = accumarray([edges(:,1); edges(:,2)], [nodes(edges(:,2),2); nodes(edges(:,1),2)], [numNodes, 1]);
                if size(nodes, 2) == 3
                    neighborSumZ = accumarray([edges(:,1); edges(:,2)], [nodes(edges(:,2),3); nodes(edges(:,1),3)], [numNodes, 1]);
                end
                
                neighborCount = accumarray([edges(:,1); edges(:,2)], 1, [numNodes, 1]);
        
                % Smooth only the interior nodes (boundary nodes stay fixed)
                for k = 1:numNodes
                    if ~boundaryNodes(k) && neighborCount(k) > 0
                        smoothedNodes(k, 1) = neighborSumX(k) / neighborCount(k);
                        smoothedNodes(k, 2) = neighborSumY(k) / neighborCount(k);
                        if size(nodes, 2) == 3
                            smoothedNodes(k, 3) = neighborSumZ(k) / neighborCount(k);
                        end
                    end
                end
            end
        end
        function boundaryNodes = detectBoundaryNodes(obj, connectivityArray, numNodes)
            % Initialize boundaryNodes as false
            boundaryNodes = false(numNodes, 1);
        
            % Get all edges from the connectivityArray
            edges = [connectivityArray(:, [1, 2]); connectivityArray(:, [2, 3]); connectivityArray(:, [3, 1])];
        
            % Sort each edge's vertex indices to ensure consistent ordering
            sortedEdges = sort(edges, 2);
            
            % Get unique edges and counts of how often each appears
            [uniqueEdges, ~, ic] = unique(sortedEdges, 'rows');
            edgeCounts = accumarray(ic, 1);
            
            % Boundary edges appear only once, get their indices
            boundaryEdgeIdx = find(edgeCounts == 1);
            
            % Mark nodes on boundary edges as boundary nodes
            boundaryEdges = uniqueEdges(boundaryEdgeIdx, :);
            boundaryNodes(unique(boundaryEdges(:))) = true;
        end
        function elems = addRectMeshTriangular2D( obj, mode, x1, y1, dx, dy, nx, ny )
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
            elems = obj.merge( newNodes, newElems );
        end
        function elems = addRectMeshTetrahedral3D( obj, mode, x, dx, nx )
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
            elems = obj.merge( newNodes, newElems );
        end
        function elems = addRectMesh3D( obj, x1, y1, z1, dx, dy, dz, nx, ny, nz, lnodes )
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
            elems = obj.merge(newNodes, newElems);
        end
        function elems = addShapedMesh2D( obj, sf, x, ndiv, pattern )
            baseMesh = Mesh();
            baseElems = baseMesh.addRectMesh2D( -1, -1, 2, 2, ndiv, ndiv, pattern );
            baseMesh.nodes = sf.computeValue( baseMesh.nodes ) * x;
            elems = obj.merge( baseMesh.nodes, baseElems );
        end
        function elems = addObjectMesh2D( obj, so, ndiv1, ndiv2, pattern )
            baseMesh = Mesh();
            baseElems = baseMesh.addRectMesh2D( -1, -1, 2, 2, ndiv1, ndiv2, pattern );
            baseMesh.nodes = so.computeValue( baseMesh.nodes );
            elems = obj.merge( baseMesh.nodes, baseElems );
        end
        function elems = addObjectMesh3D( obj, so, ndiv1, ndiv2, ndiv3, lnodes )
            baseMesh = Mesh();
            baseElems = baseMesh.addRectMesh3D( -1, -1, -1, 2, 2, 2, ndiv1, ndiv2, ndiv3, lnodes );
            baseMesh.nodes = so.computeValue( baseMesh.nodes );
            elems = obj.merge( baseMesh.nodes, baseElems );
        end
        function elems = addShapedMesh3D( obj, sf, x, ndiv, pattern )
            baseMesh = Mesh();
            baseElems = baseMesh.addRectMesh3D( -1, -1, -1, 2, 2, 2, ndiv(1), ndiv(2), ndiv(3), pattern );
            baseMesh.nodes = sf.computeValue( baseMesh.nodes ) * x;
            elems = obj.merge( baseMesh.nodes, baseElems );
        end
        function elems = addRectWithHoleInCornerMesh2D( obj, ri, x00, y00, s, div, pattern )
            xs = x00 - s;
            ys = y00 + s;
            alpha = 240*pi/180; 
            beta  = 315*pi/180;
            x1 = xs+ri*cos(alpha);
            mesh = Mesh();
            elems = mesh.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            mesh.nodes = [ 0.5*( 1-mesh.nodes(:,2) ).*((x00-x1)/2*mesh.nodes(:,1)+(x00+x1)/2)+0.5*(mesh.nodes(:,2)+1).*(xs+ri*cos( (beta-alpha)/2 .* mesh.nodes(:,1) + (beta+alpha)/2)) ...
                           0.5*( 1-mesh.nodes(:,2) ).*(y00)+0.5*(mesh.nodes(:,2)+1).*(ys+ri*sin( (beta-alpha)/2 .* mesh.nodes(:,1) + (beta+alpha)/2)) ];
            elems = [ elems; obj.merge( mesh.nodes, elems ) ];
            mesh = Mesh();
            elems1 = mesh.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            alpha1 = alpha;           
            alpha = -45*pi/180;
            beta  = pi/6;         
            y1 = ys+ri*sin(beta);
            mesh.nodes = [ 0.5*( mesh.nodes(:,1)+1 ).*x00+0.5*(1-mesh.nodes(:,1)).*(xs+ri*cos( (beta-alpha)/2 .* mesh.nodes(:,2) + (beta+alpha)/2)) ...
                           0.5*( mesh.nodes(:,1)+1 ).*((y1-y00)/2*mesh.nodes(:,2)+(y00+y1)/2)+0.5*(1-mesh.nodes(:,1)).*(ys+ri*sin( (beta-alpha)/2 .* mesh.nodes(:,2) + (beta+alpha)/2)) ];           
            % elems2 = obj.merge( mesh.nodes, elems );
            % elems3 = obj.addShapedMesh2D( ShapeFunctionL4(), [80 y00; x1 y00; 80 40; x1 ys+ri*sin(alpha1)], [div div], pattern );
            % elems4 = obj.addShapedMesh2D( ShapeFunctionL4(),  [xs+ri*cos(beta) ys+ri*sin(beta); x00 y1; x00-40 x00-80; x00 x00-80], [div div], pattern );
            elems = [ elems; obj.merge( mesh.nodes, elems1 ) ];
            elems = [ elems; ...
                      obj.addShapedMesh2D( ShapeFunctionL4(), [80 y00; x1 y00; 80 40; x1 ys+ri*sin(alpha1)], [div div], pattern ); ...
                      obj.addShapedMesh2D( ShapeFunctionL4(),  [xs+ri*cos(beta) ys+ri*sin(beta); x00 y1; x00-40 x00-80; x00 x00-80], [div div], pattern )];
        end
        function elems = addRectWithHoleMesh2D( obj, a, x0, y0, holefactor, div, pattern )
            r=a*holefactor;
            mesh = Mesh();
            elems=mesh.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            mesh.nodes = [ 0.5*( 1-mesh.nodes(:,1)).*(x0-a)+0.5*(mesh.nodes(:,1)+1).*(x0-r*cos( pi/8 .* mesh.nodes(:,2) + pi/8)) ...
                           0.5*( 1-mesh.nodes(:,1)).*(a/2.*mesh.nodes(:,2)+y0+a/2)+0.5*(mesh.nodes(:,1)+1).*(y0+r.*sin( pi/8.*mesh.nodes(:,2)+pi/8)) ];
            mesh2 = Mesh();
            elems2=mesh2.addRectMesh2D( -1, -1, 2, 2, div, div, pattern );
            mesh2.nodes  = [ 0.5*( mesh2.nodes(:,1)+1).*(x0+a)+0.5*(1-mesh2.nodes(:,1)).*(x0+r*cos( pi/8 .* mesh2.nodes(:,2) + pi/8)) ...
                         0.5*( mesh2.nodes(:,1)+1).*(a/2.*mesh2.nodes(:,2)+y0+a/2)+0.5*(1-mesh2.nodes(:,1)).*(y0+r.*sin( pi/8.*mesh2.nodes(:,2)+pi/8)) ];
            elems=[ elems;  mesh.merge( mesh2.nodes, elems2 ) ];                              
            elems=[ elems; mesh.duplicateTransformedMeshDeg2D( [x0 y0], 90, [0 0], elems ) ];
            elems=[ elems; mesh.duplicateTransformedMeshDeg2D( [x0 y0], 180, [0 0],elems ) ];            
            obj.merge( mesh.nodes, elems );
        end
        function elems = addRing2D( obj, x0, y0 , r1, r2, nr, nfi, pattern )
             mesh = Mesh();
             elems = mesh.addRectMesh2D( r1, 0, r2-r1, 2*pi, nr, nfi, pattern );
             newNodes = [ x0+mesh.nodes(:,1).*cos( mesh.nodes(:,2) ) y0+mesh.nodes(:,1).*sin( mesh.nodes(:,2) ) ];
             elems = obj.merge( newNodes, elems );
        end
        function elems = addQuarterCircle( obj, x0 , R, nr, pattern )
             mesh1 = Mesh();
             shapeFn = ShapeFunctionL4();
             elems1 = mesh1.addShapedMesh2D( shapeFn, x0+[ 0 0; R/2 0; 0 R/2; R/2/1.41 R/2/1.41 ], nr, pattern );

             sfL2 = ShapeFunctionL2();
             sl1 = ShapeObjectRectangular(sfL2,[0 R/2; R/2/1.41 R/2/1.41 ]);
             sl2 = ShapeObjectRectangular(sfL2,[R/2/1.41 R/2/1.41; R/2 0 ]);
             sc1 = CircleObject(x0,R,90,45);
             sc2 = CircleObject(x0,R,45,0);

             ms1 = MorphSpace(sl1,sc1);
             ms2 = MorphSpace(sl2,sc2);
             elems=[ elems1; ...   
                        mesh1.addObjectMesh2D( ms1, nr, nr, pattern );...
                        mesh1.addObjectMesh2D( ms2, nr, nr, pattern ) ];
             %mesh1.addShapedMesh2D( shapeFn, x0+[ R/2 0; R 0; R/2/1.41 R/2/1.41; R/1.41 R/1.41 ], nr, pattern );
             %mesh1.addShapedMesh2D( shapeFn, x0+[ 0 R/2; R/2/1.41 R/2/1.41; 0 R; R/1.41 R/1.41  ], nr, pattern );
             elems = obj.merge( mesh1.nodes, elems );
        end
        function elems = addHalfCircle( obj, x0 , R, nr, pattern )
            elems1 = obj.addQuarterCircle( x0 , R, nr, pattern );
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            elems = [ elems1; obj.addQuarterCircle( x0 , R, nr, pattern ) ];
        end
        function elems = addCircle( obj, x0 , R, nr, pattern )
            elems = obj.addQuarterCircle( x0 , R, nr, pattern );
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            elems = [elems; obj.addQuarterCircle( x0 , R, nr, pattern )];
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            elems = [elems; obj.addQuarterCircle( x0 , R, nr, pattern )];
            obj.transformMeshDeg2D( x0, 90, [0 0] );
            elems = [elems; obj.addQuarterCircle( x0 , R, nr, pattern )];
        end
        function elems = addQuarterCylinder( obj, x0 , R, h, nr, pattern )
             mesh1 = Mesh();
             shapeFn = ShapeFunctionL8();
             elems = mesh1.addShapedMesh3D( shapeFn, x0+[ 0 0 0; R/2 0 0; 0 R/2 0; R/2/1.41 R/2/1.41 0; 0 0 h; R/2 0 h; 0 R/2 h; R/2/1.41 R/2/1.41 h], [nr(1) nr(1) nr(2)], pattern );

             sfL2 = ShapeFunctionL4();
             sl1 = ShapeObjectRectangular(sfL2,x0+[0 R/2 0; R/2/1.41 R/2/1.41 0; 0 R/2 h; R/2/1.41 R/2/1.41 h ]);
             sl2 = ShapeObjectRectangular(sfL2,x0+[R/2/1.41 R/2/1.41 0; R/2 0 0; R/2/1.41 R/2/1.41 h; R/2 0 h]);
             sc1 = CylinderObject(x0,R,90,45,0,h);
             sc2 = CylinderObject(x0,R,45,0,0,h);

             ms1 = MorphSpace(sl1,sc1);
             ms2 = MorphSpace(sl2,sc2);
                
             elems = [ elems; mesh1.addObjectMesh3D( ms1, nr(1), nr(1), nr(2), pattern ) ];
             elems = [ elems; elemsmesh1.addObjectMesh3D( ms2, nr(1), nr(1), nr(2), pattern ) ];
             %mesh1.addShapedMesh2D( shapeFn, x0+[ R/2 0; R 0; R/2/1.41 R/2/1.41; R/1.41 R/1.41 ], nr, pattern );
             %mesh1.addShapedMesh2D( shapeFn, x0+[ 0 R/2; R/2/1.41 R/2/1.41; 0 R; R/1.41 R/1.41  ], nr, pattern );
             elems = obj.merge( mesh1.nodes, elems );
        end
        function elems = addHalfCylinder( obj, x0 , R, h, nr, localNodes )
            elems = obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            elems = [ elems; obj.addQuarterCylinder( x0 , R, h, nr, localNodes ) ]; 
        end
        function elems = addCylinder( obj, x0 , R, h, nr, localNodes )
            elems = obj.addQuarterCylinder( x0 , R, h, nr, localNodes );
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            elems=[elems; obj.addQuarterCylinder( x0 , R, h, nr, localNodes )];
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            elems=[elems; obj.addQuarterCylinder( x0 , R, h, nr, localNodes )];
            obj.transformMesh3DDegXY( x0, 90, [0 0 0] );
            elems=[elems; obj.addQuarterCylinder( x0 , R, h, nr, localNodes )];
        end
        function elems = addPipe3D(obj,x0,r,R,al1,al2,h1,h2,nr,nc,nz,lnodes)
            mesh=Mesh();
            elems = mesh.addRectMesh3D( r, deg2rad(al1), h1, R-r, deg2rad(al2-al1), h2-h1, nr, nc, nz, lnodes);
            mesh.transformToCylindrical3D(x0);
            elems = obj.merge(mesh.nodes,elems);
        end
        function elems = addrectPipe(obj,w,h,l1,th,nth,lnodes)
            nx=round(l1/th)*nth;
            ny=round((w-2*th)/th)*nth;
            nz=round((h-2*th)/th)*nth;
            elems = [   obj.addRectMesh3D( 0, 0, 0, l1, th, th, nx, nth, nth, lnodes );...
                        obj.addRectMesh3D( 0, 0, h-th, l1, th, th, nx, nth, nth, lnodes );...
                        obj.addRectMesh3D( 0, w-th, h-th, l1, th, th, nx, nth, nth, lnodes );...
                        obj.addRectMesh3D( 0, w-th, 0, l1, th, th, nx, nth, nth, lnodes );...
                        obj.addRectMesh3D( 0, 0, th, l1, th, h-2*th, nx, nth, nz, lnodes );...
                        obj.addRectMesh3D( 0, w-th, th, l1, th, h-2*th, nx, nth, nz, lnodes );...
                        obj.addRectMesh3D( 0, th, 0, l1, w-2*th, th, nx, ny, nth, lnodes );...
                        obj.addRectMesh3D( 0, th, h-th, l1, w-2*th, th, nx, ny, nth, lnodes ) ];
        end
        function elems = addLshape( obj, l, h, nh, pattern )
            nl = round(l/h*nh+0.5);
            elems = [   obj.addRectMesh2D( 0, 0, h, h, nh, nh, pattern ); ...
                        obj.addRectMesh2D( 0, h, h, l-h, nh, nl-nh, pattern ); ...
                        obj.addRectMesh2D( h, 0, l-h, h, nl-nh, nh, pattern )   ];
        end
        function elems = addLshape3D( obj, l, h, nh, lnodes )
            nl = round(l/h*nh+0.5);
            elems = [   obj.addRectMesh3D( 0, 0, 0, h,   h, h,   nh,    nh, nh,    lnodes);
                        obj.addRectMesh3D( 0, 0, h, h,   h, l-h, nh,    nh, nl-nh, lnodes);
                        obj.addRectMesh3D( h, 0, 0, l-h, h, h,   nl-nh, nh, nh,    lnodes) ];
        end        
        function elems = addHframe( obj, nspan, lspan, nfloor, hfloor, xp )
                [X,Y]=meshgrid(xp(1):lspan:(xp(1)+nspan*lspan),xp(2):hfloor:(xp(2)+nfloor*hfloor));
                obj.nodes = [ obj.nodes; [ X(:) Y(:) ] ];
                elems=[];
                beams=[2:nfloor+1; (2:nfloor+1)+nfloor+1]';
                cols=[1:nfloor; (1:nfloor)+1]';
                for k=1:nspan
                    elems =[ elems; (k-1)*(nfloor+1)+beams ];
                end
                for k=1:nspan    
                    elems =[ elems; (k-1)*(nfloor+1)+cols ];
                end
                elems =[ elems; nspan*(nfloor+1)+cols ];
        end
        function elems = duplicateTransformedMeshDeg2D( obj, x0, angleDeg, xm , oldelems)
            angleRad = angleDeg*pi/180;
            elems = obj.merge( [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                                 x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad)] + xm, oldelems );
                  
        end
        function elems = duplicateTransformedMeshDeg3D( obj, x0, angleDeg, xm, oldelems )
            angleRad = angleDeg*pi/180;
            elems = obj.merge( [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                         x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad)...
                         obj.nodes(:,3)] + xm, oldelems );
                  
        end
        function elems = array( obj, coord, n, oelems )
                m = max(obj.nodes)-min(obj.nodes);
                mv=zeros(1,size(obj.nodes,2));
                belems=oelems;
                bnodes=obj.nodes;
                elems=oelems;
                for k=1:n
                    mv(coord)=k*m(coord);
                    newnodes=bnodes+mv;
                    elems=[ elems; obj.merge( newnodes, belems ) ];
                end
        end
        function obj = transformMeshDeg2D( obj, x0, angleDeg, xm )
            angleRad = angleDeg*pi/180;
            obj.nodes = [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                          x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad)] + xm;
                  ;
        end
        function obj = transformMesh3DDegXY( obj, x0, angleDeg, xm )
            angleRad = angleDeg*pi/180;
            obj.nodes = [ x0(1)+(obj.nodes(:,1)-x0(1)).*cos(angleRad)-(obj.nodes(:,2)-x0(2)).*sin(angleRad)...
                          x0(2)+(obj.nodes(:,1)-x0(1)).*sin(angleRad)+(obj.nodes(:,2)-x0(2)).*cos(angleRad) obj.nodes(:,3)] + xm;
                  
        end
        function elems = transformToPolar2D( obj, x0, y0, oldelems )
             newNodes = [ x0+obj.nodes(:,1).*cos( obj.nodes(:,2) ) y0+obj.nodes(:,1).*sin( obj.nodes(:,2) ) ];
             [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
             obj.nodes = newNodes( si1, : );
             elems = si2(oldelems);
        end
        function elems = transformToCylindrical3D( obj, x0, oldelems )
             newNodes = [ x0(1)+obj.nodes(:,1).*cos( obj.nodes(:,2) ) x0(2)+obj.nodes(:,1).*sin( obj.nodes(:,2) ) obj.nodes(:,3) ];
             [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
            obj.nodes = newNodes( si1, : );
            elems = si2(oldelems);
        end
        function elems = transformToCylindricalEliptic3D( obj, x0, a, b, oldelems )
             newNodes = [ x0(1)+a*obj.nodes(:,1).*cos( obj.nodes(:,2) ) x0(2)+b*obj.nodes(:,1).*sin( obj.nodes(:,2) ) obj.nodes(:,3) ];
             [~,si1,si2] = unique( round(newNodes .* obj.tolerance), 'rows', 'stable' );
            obj.nodes = newNodes( si1, : );
            elems = si2(oldelems);
        end
        function elems = removeNodes( obj, selector, oldelems )
            nodesToRemove=find(selector.select(obj.nodes));
            nodesToLeave=find(1-selector.select(obj.nodes));
            elemsToRemove=find(sum(ismember(oldelems,nodesToRemove),2)>0);
            oldnn=size(obj.nodes,1);
            numbers=zeros(oldnn,1);
            numbers(nodesToLeave)=1:size(nodesToLeave,1);
            elems=numbers(oldelems);
            obj.nodes(nodesToRemove,:)=[];
            elems(elemsToRemove,:)=[];
        end
        function removeNodesByNumbers( obj, nodesToRemove, oldelems )
            nodesToLeave=1:size(obj.nodes,1);
            nodesToLeave(nodesToRemove)=[];
            elemsToRemove=find(sum(ismember(oldelemss,nodesToRemove),2)>0);
            oldnn=size(obj.nodes,1);
            numbers=zeros(oldnn,1);
            numbers(nodesToLeave)=1:size(nodesToLeave,1);
            elems=numbers(oldelems);
            obj.nodes(nodesToRemove,:)=[];
            elems(elemsToRemove,:)=[];
        end
        function elems = removeElemsByNumbers( obj, elemsToRemove, oldelems )
            elemsToLeave=(1:size(oldelems,1))';
            elemsToLeave(elemsToRemove)=[];
            nrem=oldelems(elemsToLeave,:);
            nodesToLeave=unique(nrem(:));
            nodesToRemove=(1:size(obj.nodes,1))';
            nodesToRemove(nodesToLeave)=[];
            oldnn=size(obj.nodes,1);
            mapNodeNumbers=zeros(oldnn,1);
            mapNodeNumbers(nodesToLeave)=1:size(nodesToLeave,1);
            elems=mapNodeNumbers(oldelems);
            obj.nodes(nodesToRemove,:)=[];
            elems(elemsToRemove,:)=[];
        end
        function transformNodesXY( obj, transformFn )
            obj.nodes = transformFn( obj.nodes );
        end
        function elems = importFEMesh( obj, mesh )
            obj.nodes=mesh.Nodes';
            elems=mesh.Elements';
        end
        function exportMeshToFile(obj, filenamebase, elems)
            save([filenamebase '_mesh.mat'],"nodes","elems");
            dlmwrite([filenamebase '_nodes.txt'],obj.nodes,'delimiter','\t','precision','%7.3f');
            dlmwrite([filenamebase '_elems.txt'],elems,'delimiter','\t','precision',6,'-append');
        end
        function plot(marker, color)
        end
    end
end

