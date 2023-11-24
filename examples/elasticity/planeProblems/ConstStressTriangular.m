clear;
close all;
res = 10;
l = 3;
tic
sf = ShapeFunctionT3;
mesh = Mesh();
mesh.addRectMeshTriangular2D( 'quad', 0, 0, 2*l, l, 2*res, res );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1))<0.001 ) );
loadedEdgeSelector = Selector( @(x)( abs(x(:,1) - 2*l)<0.001 ) );

fe=PlaneStressElem( sf, mesh.elems );
fe.plot(mesh.nodes);
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.setMaterial( material );
plot(mesh.nodes(:,1),mesh.nodes(:,2),'.');

analysis = LinearElasticity( fe, mesh );
analysis.elementLoadLineIntegral( "global", loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + [ -100 0 ] ));
analysis.fixNodes( fixedEdgeSelector, ["ux" ] );
analysis.fixClosestNode([0 0],["uy"], 0);
analysis.plotCurrentLoad();
analysis.plotSupport();

analysis.printProblemInfo();

analysis.solve();

analysis.plotMaps(["uy" "ux" "exx" "sxx" "syy" "sxy" "sHM"],0.1);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);
