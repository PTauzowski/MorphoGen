classdef CircleObject < SpatialObject
   
    
    properties
        x0,R,fi1,fi2;
    end
    
    methods
        function obj = CircleObject(x0,R,fi1,fi2)
            obj = obj@SpatialObject(20);
            obj.x0=x0;
            obj.R=R;
            obj.fi1=pi*fi1/180;
            obj.fi2=pi*fi2/180;
           
        end

        function xp = computeValue( obj, x )
            xp=zeros(size(x,1),2);
            fi = (obj.fi2-obj.fi1)/2.*x+(obj.fi1+obj.fi2)/2;
            xp(:,1) = obj.x0(1)+obj.R*cos(fi);
            xp(:,2) = obj.x0(2)+obj.R*sin(fi);
        end

        function grad = computeGradient( obj, x )
            grad=zeros(size(x,1));
            fi = (obj.fi2-obj.fi1)/2.*x+(obj.fi1+obj.fi2);
            dfi = (obj.fi2-obj.fi1)/2;
            grad(:,1) = -obj.R*sin(fi)*dfi;
            grad(:,2) = obj.R*cos(fi)*dfi;
        end

        function paint(obj)
            x=((0:obj.ndiv)/obj.ndiv*2-1)';
            y=obj.computeValue( x );
            line(y(:,1),y(:,2),'Color','green','LineWidth',1);
        end
        

    end
end

