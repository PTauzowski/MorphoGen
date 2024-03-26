classdef SpecimenModelLinear < ModelLinear
    
    properties
        const_elems;
    end

    methods
        function obj = SpecimenModelLinear(sf,a,div,E,nu,P)
            mesh1=Mesh();
            mesh1.addRectMesh2D( 0, 0, 12, 4*a, 2*div, 4*div, sf.pattern );
            mesh1.addRectWithHoleMesh2D( a, 22, 10, 4.25/10, div, sf.pattern );
            mesh1.addRectWithHoleMesh2D( a, 22, 30, 4.25/10, div, sf.pattern );
            mesh1.addRectWithHoleMesh2D( a, 57, 10, 4.25/10, div, sf.pattern );
            mesh1.addRectWithHoleMesh2D( a, 57, 30, 4.25/10, div, sf.pattern );
            mesh1.addRectMesh2D( 32, 0, 15, 4*a, 2*div, 4*div, sf.pattern );
            mesh1.duplicateTransformedMeshDeg2D( [137 0], -90, [-40 0 ] );
            mesh2=Mesh();
            mesh2.addRectMesh2D( 67, 0, 30, 4*a, 3*div, 4*div, sf.pattern );
            mesh2.duplicateTransformedMeshDeg2D( [137 0], -90, [-40 0 ] );
            mesh2.addRectMesh2D( 97, 0, 40, 40, 4*div, 4*div, sf.pattern );
            mesh2.addRectMesh2D( 67, 40, 30, 30, 3*div, 3*div, sf.pattern );
            obj.mesh=Mesh();
            obj.mesh.merge( mesh1.nodes, mesh1.elems );
            obj.const_elems = 1:size(mesh1.elems,1)';
            obj.mesh.merge( mesh2.nodes, mesh2.elems );

            meter_factor = 1000;

            obj.mesh.nodes = obj.mesh.nodes/meter_factor;
            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );      
            mX = max(obj.mesh.nodes(:,1));

            r=0.00425; %[m]
            b = 0.08; %[m]

            Fref = 0.01*P;
            qref = Fref / 4 / pi / r;
            qref2 = Fref / b;
            tolerance = 1.0E-6;
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
            fixedEdgeSelectorX = Selector( @(x)( abs(x(:,1) - 137/meter_factor) < tolerance ) );
            fixedEdgeSelectorY = Selector( @(x)( abs(x(:,2)) < tolerance) );
            loadedEdgeSelectorX = Selector( @(x)( abs(x(:,1)) < tolerance) );
            loadedEdgeSelectorY = Selector( @(x)( abs(x(:,2) - 137/meter_factor) < tolerance ) );
            holeSelector1 = Selector( @(x)( (((x(:,1) - 22/meter_factor).^2 +  (x(:,2) - 10/meter_factor).^2 )-(r)^2 ) < tolerance ) );
            holeSelector2 = Selector( @(x)( (((x(:,1) - 22/meter_factor).^2 +  (x(:,2) - 30/meter_factor).^2 )-(r)^2 )  < tolerance  ) );
            holeSelector3 = Selector( @(x)( (((x(:,1) - 57/meter_factor).^2 +  (x(:,2) - 10/meter_factor).^2 )-(r)^2 )  < tolerance ) );
            holeSelector4 = Selector( @(x)( (((x(:,1) - 57/meter_factor).^2 +  (x(:,2) - 30/meter_factor).^2 )-(r)^2 )  < tolerance ) );
            holeSelector5 = Selector( @(x)( (((x(:,1) - 107/meter_factor).^2 + (x(:,2) - 115/meter_factor).^2 )-(r)^2 )  < tolerance ) );
            holeSelector6 = Selector( @(x)( (((x(:,1) - 127/meter_factor).^2 + (x(:,2) - 115/meter_factor).^2 )-(r)^2 )  < tolerance ) );
            holeSelector7 = Selector( @(x)( (((x(:,1) - 107/meter_factor).^2 + (x(:,2) - 80/meter_factor).^2 )-(r)^2 )  < tolerance ) );
            holeSelector8 = Selector( @(x)( (((x(:,1) - 127/meter_factor).^2 + (x(:,2) - 80/meter_factor).^2 )-(r)^2 )  < tolerance ) );
            
            %obj.analysis.plotSelectedNodeNumbers( loadedEdgeSelector );
            
            obj.analysis.elementLoadLineIntegral( "global", loadedEdgeSelectorX, "ux", @(x)( x(:,1)*0 - qref2 ));
            obj.analysis.elementLoadLineIntegral( "global", loadedEdgeSelectorY, "uy", @(x)( x(:,2)*0 + qref2 ));
            
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector1, "ux", @(x)( x(:,1)*0 - qref ));
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector2, "ux", @(x)( x(:,1)*0 - qref ));
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector3, "ux", @(x)( x(:,1)*0 - qref ));
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector4, "ux", @(x)( x(:,1)*0 - qref ));
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector5, "uy", @(x)( x(:,2)*0 + qref ));
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector6, "uy", @(x)( x(:,2)*0 + qref ));
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector7, "uy", @(x)( x(:,2)*0 + qref ));
            % obj.analysis.elementLoadLineIntegral( "global", holeSelector8, "uy", @(x)( x(:,2)*0 + qref ));
            obj.analysis.fixNodes( fixedEdgeSelectorX, "ux"); 
            obj.analysis.fixNodes( fixedEdgeSelectorY, "uy" );

            obj.result_number=17;
            
        end

        function setupLoad(obj,P)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(P(1), 0.3);
            obj.fe.setMaterial( material );   
        end
       
    end
end
