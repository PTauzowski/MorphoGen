clear;
close all;
res = 20;
l = 3;
tic
sfL4 = ShapeFunctionL4;
mesh = Mesh();
elems=mesh.addRectMesh2D(0, 0, 2*l, l, 2*res, res, sfL4.pattern);
fe=PlaneStressElem( sfL4, elems );
mX = max(mesh.nodes(:,1));
fe.props.h=1;

material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.setMaterial( material );

analysis = LinearElasticity( fe, mesh );

fixedEdgeSelector = Selector( @(x)( abs(x(:,1))<0.001 ) );
loadedEdgeSelector = Selector( @(x)( abs(x(:,1) - 2*l)<0.001 ) );
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysis.elementLoadLineIntegral( "global", loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + [75 -75] ));
analysis.printProblemInfo();

mesh.transformMeshDeg2D( [0 0], 45, [0 0] );
fe.plot(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();

q = analysis.solve();

analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.2);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);

