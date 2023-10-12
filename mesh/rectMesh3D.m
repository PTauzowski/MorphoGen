function [ nodes, elems ] = rectMesh3D( nx, ny, nz, x0, y0, z0, dx, dy, dz, pattern )

    dim = max(max( pattern ));
    ne  = size( pattern, 1 ); % number of nodes in single element
   
    x = linspace(x0,x0+dx,nx+1);
    y = linspace(y0,y0+dy,ny+1);
    z = linspace(z0,z0+dz,nz+1);
    
    [X,Y,Z] = meshgrid(x,y,z);
    
      
    nodes   = [ X(:) Y(:) Z(:) ];
    nn      = size(nodes,1);
    nfe     = nx * ny * nz;
    
    elems = zeros( nfe, ne );
    grid  = X;
    grid(:) = 1:nn;
    count = 1;
    for k=1:ny
      for l=1:nx
         for m=1:nz
             for n = 1:ne
                elems(count,n) =  grid( k*dim + pattern(n,1),  l*dim + pattern(n,2), m*dim + pattern(n,3) );
             end   
             count = count + 1;
         end
      end
    end
      
end