classdef CylinderObject < SpatialObject
   
    
    properties
        x0,R,fi1,fi2,z1,z2;
    end
    
    methods
        function obj = CylinderObject(x0,R,fi1,fi2,z1,z2)
            obj = obj@SpatialObject(20);
            obj.x0=x0;
            obj.R=R;
            obj.fi1=pi*fi1/180;
            obj.fi2=pi*fi2/180;
            obj.z1=z1;
            obj.z2=z2;
           
        end

        function xp = computeValue( obj, x )
            xp=zeros(size(x,1),3);
            fi = (obj.fi2-obj.fi1)/2.*x(:,1)+(obj.fi1+obj.fi2)/2;
            xp(:,1) = obj.x0(1)+obj.R*cos(fi);
            xp(:,2) = obj.x0(2)+obj.R*sin(fi);
            xp(:,3) = (obj.z2-obj.z1)/2.*x(:,2)+(obj.z1+obj.z2)/2;
        end

        function grad = computeGradient( obj, x )
            grad=zeros(size(x,1),3,2);
            fi = (obj.fi2-obj.fi1)/2.*x(:,1)+(obj.fi1+obj.fi2);
            dfi = (obj.fi2-obj.fi1)/2;
            grad(:,1,1) = -obj.R*sin(fi)*dfi;
            grad(:,2,1) = obj.R*cos(fi)*dfi;
            grad(:,3,1) = 0;
            grad(:,1,2) = 0;
            grad(:,2,2) = 0;
            grad(:,3,2) = (obj.z2-obj.z1)/2;
        end

        function paint(obj)
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
            %patch('Vertices', ye, 'Faces', 1:size(ye,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','None');
            if  size(y1,2)==2
                line(reshape(y1(:,1),obj.ndiv+1,obj.ndiv+1),reshape(y1(:,2),obj.ndiv+1,obj.ndiv+1),'Color','blue');
                line(reshape(y2(:,1),obj.ndiv+1,obj.ndiv+1),reshape(y2(:,2),obj.ndiv+1,obj.ndiv+1),'Color','blue');
            elseif size(y1,2)==3
                line(reshape(y1(:,1),obj.ndiv+1,obj.ndiv+1),reshape(y1(:,2),obj.ndiv+1,obj.ndiv+1),reshape(y1(:,3),obj.ndiv+1,obj.ndiv+1),'Color','blue');
                line(reshape(y2(:,1),obj.ndiv+1,obj.ndiv+1),reshape(y2(:,2),obj.ndiv+1,obj.ndiv+1),reshape(y2(:,3),obj.ndiv+1,obj.ndiv+1),'Color','blue');
            end
        end
        

    end
end

