classdef MorphSpace < Function
    
    properties
        ss1,ss2,ndiv
    end
    
    methods
        function obj = MorphSpace(ss1,ss2)
            obj=obj@Function(ss1.dim+1,0.0001)
            obj.ss1 = ss1;
            obj.ss2 = ss2;
            obj.ndiv=10;
        end
        
        function xp = computeValue( obj, x )
            xp=(1-x(:,1))/2.*obj.ss1.computeValue(x(:,2:end))+(1+x(:,1))/2.*obj.ss2.computeValue(x(:,2:end));
        end

        function grad = computeGradient( obj, x )
            grad=zeros(size(x,1),size(x,2),size(x,2));
            grad(:,:,1)=(obj.ss2.computeValue(x(:,2:end))-obj.ss1.computeValue(x(:,2:end)))/2;
            grad(:,:,2:end)=(1-x(:,1))/2.*obj.ss1.computeGradient(x(:,2:end))+(1+x(:,1))/2.*obj.ss2.computeGradient(x(:,2:end));
            
        end

        function paint(obj)
            switch size(obj.ss1.shapeFunctions.localNodes,2)
                case 1
                   x=((0:obj.ndiv)/obj.ndiv*2-1)';
                   xy1=obj.ss1.computeValue( x );
                   xy2=obj.ss2.computeValue( x );
                   obj.ss1.paint();
                   obj.ss2.paint();
                   line([xy1(:,1) xy2(:,1)]',[xy1(:,2) xy2(:,2)]','Color','blue','LineWidth',0.25);
                case 2
                    [Y,X]=meshgrid(((0:obj.ndiv)/obj.ndiv*2-1),((0:obj.ndiv)/obj.ndiv*2-1));
                    x=[X(:) Y(:)];
                    is_in_border = [ find(any(x == -1,2)); find(any(x == 1,2))];
                    xy1=obj.ss1.computeValue( x(is_in_border,:) );
                    xy2=obj.ss2.computeValue( x(is_in_border,:) );
                    obj.ss1.paint();
                    obj.ss2.paint();
                    if size(xy1,2)==2
                        line([xy1(:,1) xy2(:,1)]',[xy1(:,2) xy2(:,2)]','Color','blue','LineWidth',0.25);
                    elseif size(xy1,2)==3
                        line([xy1(:,1) xy2(:,1)]',[xy1(:,2) xy2(:,2)]', [xy1(:,3) xy2(:,3)]','Color','blue','LineWidth',0.25);
                    end
                case 3
                    disp('3D');
            otherwise
                disp('wrong ');
            end
        end
    end
end

