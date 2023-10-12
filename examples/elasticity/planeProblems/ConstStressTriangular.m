clear;
close all;
res = 10;
l = 3;
tic
sfL9 = ShapeFunctionT3;
mesh = Mesh();
mesh.addRectMeshTriangular2D( 'quad', 0, 0, 2*l, l, 2*res, res );
fe=PlaneStressElem( sfL9, mesh.elems );
fe.plotSolid(mesh.nodes);
fe.props.h=1;
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.setMaterial( material );
plot(mesh.nodes(:,1),mesh.nodes(:,2),'.');

analysis = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( @(x)( x(:,1) ) );
loadedEdgeSelector = Selector( @(x)( x(:,1) - 2*l ) );
circleSelector = Selector( @(x)( ((x(:,1) - 1.5).^2 + (x(:,2) - 1.5).^2 )-1.5 ), 0.2 );
%fedges = fe.findEdges( problem.findNodes( loadedEdgeSelector ));
%problem.plotSelectedNodeNumbers( loadedEdgeSelector );
%P = problem.loadEdgesGlobal( loadedEdgeSelector, "ux", @(x)( x*0 + [-150 0 ] ));
analysis.elementLoadLineIntegral( "global", loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + [ -100 0 ] ));
analysis.fixNodes( fixedEdgeSelector, ["ux" ] );
analysis.fixClosestNode([0 0],["uy"], 0);
analysis.plotCurrentLoad();
analysis.plotSupport();

analysis.printProblemInfo();

analysis.solve();

analysis.plotMaps(["uy" "ux" "exx" "sxx" "syy" "sxy" "sHM"],0.1);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);
