classdef ShapeSpace < Function
   
    properties
        shapeFunctions;
        points;
        ndiv;
    end
    
    methods
        
        function obj = ShapeSpace(shapeFunctions,points)
            obj.shapeFunctions = shapeFunctions;
            obj.points=points;
            obj.ndiv=10;
        end
        
        function xp = computeValue( obj, x )
            N=obj.shapeFunctions.computeValue(x);
            xp=zeros(size(x,1),size(x,2));
            for k=1:size(obj.points,2)
                xp(:,k)=N*obj.points(:,k);
            end
        end

        function grad = computeGradient( obj, x )
            dN=obj.shapeFunctions.computeValue(x);
            grad=zeros(size(x,1),size(x,2),size(x,2));
            for d=1:size(x,2)
                for k=1:size(x,2)
                    grad(:,k,d)=dN(:,:,d)*obj.points(:,k);
                end
            end
        end

        function plot(obj)
            switch size(obj.shapeFunctions.localNodes,2)
                case 1
                   x=((0:obj.ndiv)/obj.ndiv*2-1)';
                   y=computeValue( obj, x );
                   plot(y(:,1),y(:,2));
                case 2
                    [Y,X]=meshgrid(((0:obj.ndiv)/obj.ndiv*2-1),((0:obj.ndiv)/obj.ndiv*2-1));
                    x1=[X(:) Y(:)];
                    x2=[Y(:) X(:)];
                    y1=computeValue( obj, x1 );
                    y2=computeValue( obj, x2 );
                    eb=[((0:obj.ndiv)/obj.ndiv*2-1)' ((0:obj.ndiv)/obj.ndiv*2-1)'];
                    er=[((0:obj.ndiv)/obj.ndiv*2-1)' ((0:obj.ndiv)/obj.ndiv*2-1)'];
                    eu=[((0:obj.ndiv)/obj.ndiv*2-1)' ((0:obj.ndiv)/obj.ndiv*2-1)'];
                    el=[((0:obj.ndiv)/obj.ndiv*2-1)' ((0:obj.ndiv)/obj.ndiv*2-1)'];
                    eb(:,2)=-1;
                    er(:,1)=1;
                    eu(:,2)=1;
                    el(:,1)=-1;
                    xe=[eb; er; eu(obj.ndiv+1:-1:1,:); el(obj.ndiv+1:-1:1,:)];
                    ye=obj.computeValue( xe );
                    patch('Vertices', ye, 'Faces', 1:size(ye,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','None');
                    line(reshape(y1(:,1),obj.ndiv+1,obj.ndiv+1),reshape(y1(:,2),obj.ndiv+1,obj.ndiv+1),'Color','blue');
                    line(reshape(y2(:,1),obj.ndiv+1,obj.ndiv+1),reshape(y2(:,2),obj.ndiv+1,obj.ndiv+1),'Color','blue');
                case 3
                    disp('3D');
            otherwise
                disp('wrong ');
            end
        end
    end
end

