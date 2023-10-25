clear;
close all;
res = 20;
l = 3;
conn=l/res*10;
tic
sf = ShapeFunctionL16;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, 2*l, l, 2*res, res, sf.pattern);
mesh2 = Mesh();
mesh2.addRectMesh2D(0, l, 2*l, l, 2*res, res, sf.pattern);
elems1 = mesh.elems;
connectionSelector = Selector( @(x)( ( abs( x(:,2) - l ) < 1.0E-04) & (( x(:,1) <= conn ) | ( 2*l - x(:,1) <= conn )  )  ) );
elems2 = mesh.connect( connectionSelector, mesh2.nodes, mesh2.elems );

fe1=PlaneStressElem( sf, elems1 );
fe1.plotSolid(mesh.nodes);
fe1.props.h=1;
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe1.setMaterial( material );

fe2=PlaneStressElem( sf, elems2 );
fe2.plotSolid(mesh.nodes);
fe2.props.h=1;
material = PlaneStressMaterial('mat2');
material.setElasticIzo(81000, 0.3);
fe2.setMaterial( material );

analysis = LinearElasticity( {fe1 fe2}, mesh );
loadedEdgeSelector1 = Selector( @(x)( abs(x(:,2))<0.001 ) );
loadedEdgeSelector2 = Selector( @(x)( abs(x(:,2) - 2*l ) <0.001 ) );
analysis.elementLoadLineIntegral( "global", loadedEdgeSelector1, ["ux" "uy"], @(x)( x*0 + [0  -100] ));
analysis.elementLoadLineIntegral( "global", loadedEdgeSelector2, ["ux" "uy"], @(x)( x*0 + [0   100] ));

analysis.fixClosestNode( [ 0 l], ["ux" "uy"], [0 0] );
analysis.fixClosestNode( [ 2*l l ],["uy"], [0] );

analysis.printProblemInfo();
fe1.plotSolid(mesh.nodes);
fe2.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();
analysis.plotSelectedNodes( connectionSelector );

q = analysis.solve();

analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.03);
fe1.plotWired(mesh.nodes,analysis.qnodal,0.05);
fe2.plotWired(mesh.nodes,analysis.qnodal,0.05);
