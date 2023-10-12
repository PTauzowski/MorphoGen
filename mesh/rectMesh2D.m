function [ nodes, elems ] = rectMesh2D( nx, ny, x0, y0, dx, dy, pattern )

    dim = max(max( pattern ));
    nn = ( dim *  nx + 1 ) * ( dim *  ny + 1 );
    ne = size( pattern, 1 ); % number of nodes in single element
    nfe = nx * ny;
    nodes = zeros( nn, 2 );
    elems = zeros( nfe, ne );
    nodes(1:nn,1) = x0 + dx * rem( (1:nn) - 1 , nx + 1 ) / nx;
    nodes(1:nn,2) = y0 + dy * ( floor( ( (1:nn) - 1 ) / ( nx + 1 ) ) / ny );
    gr = 1:nn;
    grid = reshape( gr, dim *  nx + 1 , dim * ny + 1  );
      
    for k=1:nfe
      ix = dim * rem((k)-1,nx)+1;
      iy = dim * floor(((k)-1)/(nx))+1;
      for i=1:ne
        elems(k,i) = grid( ix + pattern(i,1),  iy + pattern(i,2) ) ;
      end
    end
      
end