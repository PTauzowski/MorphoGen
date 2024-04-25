classdef Pylon2DModel < ModelLinearLoad
    
    methods
        function obj = Pylon2DModel(sf,h,h1,b,b1,resb,E,nu,xp,resX)
            resh=resb*round(h/b);
            resh1=resb*round(h/h1);
            resb1=resb*round(b1/b);

            x1=-b1/2;
            x2=-b/2;
            x3=b/2;
            x4=b1/2;
            z1=0;
            z2=h-h1;
            z3=h;


            obj.mesh = Mesh();
            obj.mesh.addRectMesh2D(x2, 0, b,  h-h1, resb, resh-resh1, sf.pattern );
            obj.mesh.addRectMesh2D(x2, z2, b, h1, resb, resh1, sf.pattern );
            obj.mesh.addRectMesh2D(x1, z2, x2-x1, h1, round((resb1-resb)/2), resh1, sf.pattern );
            obj.mesh.addRectMesh2D(x3, z2, x4-x3, h1, round((resb1-resb)/2), resh1, sf.pattern );
            

            fixedFaceSelector = Selector( @(x)( abs(x(:,2) - l)<0.001 ) );
            loadedFaceSelector = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );
            
            obj.fe = PlaneStressElem( sf, obj.mesh.elems );
            
            obj.fe.plot(obj.mesh.nodes);
            obj.fe.props.h=1;
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material)
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            
            obj.xp=xp;
            % obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], [1 0 0]);
            % obj.analysis.createNextRightHandSideVector();
            % obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], [0 1 0]);
            % obj.analysis.createNextRightHandSideVector();
            % obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], [0 0 1]);
            % obj.analysis.createNextRightHandSideVector();

            obj.analysis.loadClosestNode([-b1/2   h-h1],["ux" "uy" ], [0 -1]);
            obj.analysis.loadClosestNode([ b1/2   h-h1],["ux" "uy" ], [0 -1]);
           
          
            obj.analysis.fixClosestNode([-b/2   0],["ux" "uy" ],[ 0 0]);
            obj.analysis.fixClosestNode([ b/2   0],["ux" "uy" ],[ 0 0]);
            
            obj.analysis.printProblemInfo();

            obj.setResultNode(resX);
            obj.result_number=25;
            obj.setX(ones(obj.analysis.getTotalElemsNumber(),1));
        end

         function setupVariables(obj,E,nu,P)
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], P);
         end

         function pstress = computePenalizedHMstress(obj,E,nu,pressure,penalty)
            obj.analysis.clearCurrentLoad();
            obj.setupVariables(E,nu,pressure);
            obj.analysis.solveWeighted(obj.x);
            obj.analysis.computeElementResults(obj.x);
            pstress=sum(obj.fe.results.nodal.all(:,14).*obj.fe.results.nodal.all(:,obj.result_number))^(1/penalty);
            %pstress=sum(obj.fe.results.nodal.all(:,obj.result_number))^(1/penalty);
        end
       
    end
end

