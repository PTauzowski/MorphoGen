clear;
close all;
res = 20;
l = 3;
tic
sfL4 = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, 2*l, l, 2*res, res, sfL4.pattern);
fe=PlaneStressElem( sfL4, mesh.elems );
mX = max(mesh.nodes(:,1));
fe.props.h=1;

material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.setMaterial( material );

analysis = LinearElasticity( fe, mesh);

fixedEdgeSelector = Selector( @(x)( x(:,1) ) );
loadedEdgeSelector = Selector( @(x)( x(:,1) - 2*l ) );
circleSelector = Selector( @(x)( ((x(:,1) - 1.5).^2 + (x(:,2) - 1.5).^2 )-1.5 ), 0.2 );
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
%fedges = fe.findEdges( problem.findNodes( loadedEdgeSelector ));
%problem.plotSelectedNodeNumbers( loadedEdgeSelector );

analysis.elementLoadLineIntegral( "global", loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + [75 -75] ));

analysis.printProblemInfo();

mesh.transformMeshDeg2D( [0 0], 45, [0 0] );
fe.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();

q = analysis.solve();

analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.2);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);
