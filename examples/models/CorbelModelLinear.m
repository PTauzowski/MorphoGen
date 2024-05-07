classdef CorbelModelLinear < ModelLinearLoad

    properties
        fl,h,b,hb,hc,lc;
    end
    
    methods
        function obj = CorbelModelLinear(sf,b,h,lc,fl,hc,resb,E,nu)
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
            obj.mesh.addRectMesh2D(  0, 0, b, hb, resb, reshb, sf.pattern );
            obj.mesh.addRectMesh2D(  0, hb+hc, b, hb, resb, reshb, sf.pattern );
            obj.mesh.addRectMesh2D(  0, hb, b, hc, resb, reshc, sf.pattern );
            obj.mesh.addRectMesh2D(  b, hb, lc, hc, reslc, reshc, sf.pattern );
            obj.result_node = obj.mesh.findClosestNode([b+lc h/2]);

            fixedEdgeSelector1 = Selector( @(x)( abs(x(:,2)-h) ) < 0.001 );
            fixedEdgeSelector2 = Selector( @(x)( abs(x(:,2)) ) < 0.001 );

            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );            
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
            obj.setOneX();
            obj.P0fem=obj.analysis.Pfem;
        end

        function setupLoad(obj,P)
            obj.analysis.clearCurrentLoad();
            obj.analysis.loadClosestNode([obj.b+obj.fl obj.hb+obj.hc],["ux" "uy"], [0 P(1)]);
            obj.analysis.loadClosestNode([obj.b+obj.lc obj.h/2],["ux" "uy"], [P(2) 0]);
        end
       
    end
end

