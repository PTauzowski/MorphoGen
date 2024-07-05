classdef CorbelModelMultiMat < FEModel

    properties
        fl,h,b,hb,hc,lc, analysis, fe, result_node, result_number;
        mat1,mat2,mat3,nu;
    end
    
    methods
        function obj = CorbelModelMultiMat(sf,b,h,lc,fl,hc,resb,E1,E2,E3,nu)
            obj.h=h;
            obj.b=b;
            obj.hc=hc;
            obj.lc=lc;
            obj.fl=fl;
            hb=(h-hc)/2;
            obj.hb=hb;
            obj.nu=nu;
            reshb=round(resb*hb/b);
            reshc=round(resb*hc/b);
            reslc=round(resb*lc/b);
            
            obj.mesh = Mesh();    
            %addRectMesh2D( obj, x1, y1, dx, dy, nx, ny, pattern )
            elems1 = obj.mesh.addRectMesh2D(  0, 0, b, hb, resb, reshb, sf.pattern );
            elems2 = obj.mesh.addRectMesh2D(  0, hb+hc, b, hb, resb, reshb, sf.pattern );
            elems2 = [ elems2; obj.mesh.addRectMesh2D(  0, hb, b, hc, resb, reshc, sf.pattern ) ];
            elems3 = obj.mesh.addRectMesh2D(  b, hb, lc, hc, reslc, reshc, sf.pattern );
            %obj.result_node = obj.mesh.findClosestNode([b h-hb]);
            obj.result_node = obj.mesh.findClosestNode([b+lc h/2]);

            fixedEdgeSelector1 = Selector( @(x)( abs(x(:,2)-h) ) < 0.001 );
            fixedEdgeSelector2 = Selector( @(x)( abs(x(:,2)) ) < 0.001 );

            obj.fe={ PlaneStressElem( sf, elems1 ) PlaneStressElem( sf, elems2 ) PlaneStressElem( sf, elems3 ) };
            obj.mat1 = PlaneStressMaterial('mat1');
            obj.mat1.setElasticIzo(E1, nu);
            obj.mat2 = PlaneStressMaterial('mat2');
            obj.mat2.setElasticIzo(E2, nu);
            obj.mat3 = PlaneStressMaterial('mat3');
            obj.mat3.setElasticIzo(E3, nu);
            obj.fe{1}.setMaterial( obj.mat1 );   
            obj.fe{2}.setMaterial( obj.mat2 ); 
            obj.fe{3}.setMaterial( obj.mat3 ); 
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            obj.analysis.loadClosestNode([b+fl hb+hc],["ux" "uy"], [0 1]);
            obj.analysis.loadClosestNode([b+lc h/2],["ux" "uy"], [1 0]);
            obj.analysis.fixNodes( fixedEdgeSelector1, ["ux" "uy"] ); 
            obj.analysis.fixNodes( fixedEdgeSelector2, ["ux" "uy"] ); 

            obj.analysis.printProblemInfo();
            obj.analysis.initializeResults();
            obj.x=ones(obj.analysis.getTotalElemsNumber(),1);
            obj.result_number=17;
 %          obj.setOneX();
  %          obj.P0fem=obj.analysis.Pfem;
        end

        function obj = setX(obj,xopt)
            obj.x=xopt;
        end

        function setupLoad(obj,P)
            obj.analysis.clearCurrentLoad();
            obj.analysis.loadClosestNode([obj.b+obj.fl obj.hb+obj.hc],["ux" "uy"], [0 P(1)]);
            obj.analysis.loadClosestNode([obj.b+obj.lc obj.h/2],["ux" "uy"], [P(2) 0]);
        end

        function plotModelMesh( obj )
            obj.analysis.plotFiniteElements();
        end

         function obj = plotModel( obj )
            obj.analysis.plotFiniteElements();
            obj.analysis.plotCurrentLoad();
            obj.analysis.plotSupport();
            % plot(obj.mesh.nodes(obj.result_node,1),obj.mesh.nodes(obj.result_node,2),"Marker","o");
         end

          function solveWeighted(obj)
            obj.analysis.qfem=obj.analysis.solveWeighted(obj.x);
            obj.analysis.initializeResults();
            obj.analysis.computeElementResults();
          end

          function sp = computePenalizedStress(obj,penalty,points)
            np=size(points,1);
            sp=zeros(np,1);
            for k=1:np
                E1=points(k,1);
                E2=points(k,2);
                E3=points(k,3);
                P1=points(k,4);
                P2=points(k,5);
                
                obj.mat1.setElasticIzo(E1, obj.nu);
                obj.mat2.setElasticIzo(E2, obj.nu);
                obj.mat3.setElasticIzo(E3, obj.nu);

                obj.analysis.loadClosestNode([obj.b+obj.fl obj.hb+obj.hc],["ux" "uy"], [0 P1]);
                obj.analysis.loadClosestNode([obj.b+obj.lc obj.h/2],["ux" "uy"], [P2 0]);
                
                obj.analysis.computeElementResults(obj.x);
                sp(k)=sum(obj.fe{1}.results.nodal.all(:,18).*obj.fe{1}.results.nodal.all(:,obj.result_number))^(1/penalty);
            end
        end   
       
    end
end

