classdef MeshRect2D < Mesh
    
    properties
        x1, y1, dx, dy, nx, ny;
    end
    
    methods
        function obj = MeshRect2D(x1, y1, dx, dy, nx, ny, sf)
            obj = obj@Mesh(sf);
            obj.x1 = x1;
            obj.y1 = y1;
            obj.dx = dx;
            obj.dy = dy;
            obj.nx = nx;
            obj.ny = ny;
        end
        function obj = init(obj, x1, y1, dx, dy, nx, ny, sf)
            obj.sf = sf;
            obj.x1 = x1;
            obj.y1 = y1;
            obj.dx = dx;
            obj.dy = dy;
            obj.nx = nx;
            obj.ny = ny;
        end
        function generateRectangular(obj)
            dim = max(max( obj.sf.pattern ));
            nn = ( dim *  obj.nx + 1 ) * ( dim *  obj.ny + 1 );
            ne = size( obj.sf.pattern, 1 ); 
            nfe = obj.nx * obj.ny;
            obj.nodes = zeros( nn, 2 );
            obj.elems = zeros( nfe, ne );
            obj.nodes(1:nn,1) = obj.x1 + obj.dx * rem( (1:nn) - 1 , obj.nx + 1 ) / obj.nx;
            obj.nodes(1:nn,2) = obj.y1 + obj.dy * ( floor( ( (1:nn) - 1 ) / ( obj.nx + 1 ) ) / obj.ny );
            gr = 1:nn;
            grid = reshape( gr, dim *  obj.nx + 1 , dim * obj.ny + 1  );
            for k=1:nfe
              ix = dim * rem((k)-1,obj.nx)+1;
              iy = dim * floor(((k)-1)/(obj.nx))+1;
              for i=1:ne
                obj.elems(k,i) = grid( ix + obj.sf.pattern(i,1),  iy + obj.sf.pattern(i,2) ) ;
              end
            end
        end
    end
end

