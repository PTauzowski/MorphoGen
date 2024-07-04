classdef CorbelModelMultiMat < FEModel

    properties
        fl,h,b,hb,hc,lc, analysis, fe, result_node, result_number;
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
            material1 = PlaneStressMaterial('mat1');
            material1.setElasticIzo(E1, nu);
            material2 = PlaneStressMaterial('mat2');
            material2.setElasticIzo(E2, nu);
            material3 = PlaneStressMaterial('mat3');
            material3.setElasticIzo(E3, nu);
            obj.fe{1}.setMaterial( material1 );   
            obj.fe{2}.setMaterial( material2 ); 
            obj.fe{3}.setMaterial( material3 ); 
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            obj.analysis.loadClosestNode([b+fl hb+hc],["ux" "uy"], [0 1]);
            obj.analysis.createNextRightHandSideVector();
            obj.analysis.loadClosestNode([b+lc h/2],["ux" "uy"], [1 0]);
            obj.analysis.createNextRightHandSideVector();
            obj.analysis.fixNodes( fixedEdgeSelector1, ["ux" "uy"] ); 
            obj.analysis.fixNodes( fixedEdgeSelector2, ["ux" "uy"] ); 

            obj.analysis.printProblemInfo();
            obj.analysis.initializeResults();
            obj.x=ones(obj.analysis.getTotalElemsNumber(),1);
            obj.result_number=17;
 %         obj.setOneX();
  %          obj.P0fem=obj.analysis.Pfem;
        end

        function setupLoad(obj,P)
            obj.analysis.clearCurrentLoad();
            obj.analysis.loadClosestNode([obj.b+obj.fl obj.hb+obj.hc],["ux" "uy"], [0 P(1)]);
            obj.analysis.loadClosestNode([obj.b+obj.lc obj.h/2],["ux" "uy"], [P(2) 0]);
        end

        function plotModelMesh( obj )
            obj.analysis.plotFiniteElements();
        end

         function plotModel( obj )
            obj.analysis.plotFiniteElements();
            obj.analysis.plotCurrentLoad();
            obj.analysis.plotSupport();
            % plot(obj.mesh.nodes(obj.result_node,1),obj.mesh.nodes(obj.result_node,2),"Marker","o");
        end
       
    end
end

