clear;
close all;
l = 3;
x0 = 0;
y0 = 0;
r1=5;
r2=10;
div=20;
sf = ShapeFunctionL4;
mesh = Mesh();
mesh.addRing2D( x0, y0 , r1, r2, div, round(3*pi*div), sf.pattern );
%mesh.transformMeshDeg2D( [137 0], -90, [-137 0] );
fe=PlaneStressElem( sf, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.0);
fe.setMaterial( material );

fe.props.h=1;
fe.plotSolid(mesh.nodes);
plot(mesh.nodes(:,1),mesh.nodes(:,2),'.');

analysis = LinearElasticity( fe, mesh );
circleSelectorInt = Selector( @(x)( ((x(:,1) - x0).^2 + (x(:,2) - y0).^2 ) - (r1)^2 ), 0.01 );
circleSelectorOut = Selector( @(x)( ((x(:,1) - x0).^2 + (x(:,2) - y0).^2 ) - (r2)^2 ), 0.01 );
analysis.elementLoadLineIntegral( "local", circleSelectorOut,  ["ux" "uy"], @(x)( x*0 + [ 0 -15 ] ));
analysis.fixClosestNode( [ x0-r1 y0], ["ux" "uy"], [0 0] );
analysis.fixClosestNode([ x0+r1 y0],"uy", 0 );
analysis.plotCurrentLoad();
analysis.plotSupport();

analysis.printProblemInfo();

analysis.solve();
figure
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);
analysis.plotMaps(["uy" "ux" "sxx" "syy" "s1" "s2" "sxy" "sHM"],0.1);
