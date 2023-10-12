classdef MeshRect3D < Mesh
   properties
        x1, y1, z1, dx, dy, dz, nx, ny, nz;
    end
    
    methods
        function obj = MeshRect3D(x1, y1, z1, dx, dy, dz, nx, ny, nz, sf)
            obj = obj@Mesh(sf);
            obj.x1 = x1;
            obj.y1 = y1;
            obj.z1 = z1;
            obj.dx = dx;
            obj.dy = dy;
            obj.dz = dz;
            obj.nx = nx;
            obj.ny = ny;
            obj.nz = nz;
        end
        function obj = init(x1, y1, z1, dx, dy, dz, nx, ny, nz, sf)
            obj.x1 = x1;
            obj.y1 = y1;
            obj.z1 = z1;
            obj.dx = dx;
            obj.dy = dy;
            obj.dz = dz;
            obj.nx = nx;
            obj.ny = ny;
            obj.nz = nz;
            obj.sf = sf;
        end
        function obj = generateRectangular(obj)
            ne  = size( obj.sf.pattern, 1 ); 
            x = linspace(obj.x1,obj.x1+obj.dx,obj.nx+1);
            y = linspace(obj.y1,obj.y1+obj.dy,obj.ny+1);
            z = linspace(obj.z1,obj.z1+obj.dz,obj.nz+1);

            [X,Y,Z] = meshgrid(x,y,z);

            obj.nodes = [ X(:) Y(:) Z(:) ];
            nn = size(nodes,1);
            nfe = obj.nx * obj.ny * obj.nz;

            obj.elems = zeros(nfe, ne);
            grid = X;
            grid(:) = 1:nn;
            count = 1;
            for k=1:obj.ny
              for l=1:obj.nx
                 for m=1:obj.nz
                     for n = 1:ne
                        obj.elems(count,n) =  grid( k + obj.sf.pattern(n,1),  l + obj.sf.pattern(n,2), m + obj.sf.pattern(n,3) );
                     end   
                     count = count + 1;
                 end
              end
            end
        end
    end
end

