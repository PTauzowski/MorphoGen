clear;
close all;
res = 10;
l = 3;
tic
sfT3 = ShapeFunctionT3;
mesh = Mesh();
P = [0 0; 1 0; 1 0.7; 2 2; 2 3; 1 3; 1 2; 0 1];
C = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1];
mesh.addDelaunayMesh2D( P, C, 500 );
fe=PlaneStressElem( sfT3, mesh.elems );
fe.plot(mesh.nodes);
fe.props.h=1;
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.setMaterial( material );
plot(mesh.nodes(:,1),mesh.nodes(:,2),'.');

analysis = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1))<0.001 ) );
loadedEdgeSelector = Selector( @(x)( abs(x(:,1) - 2*l)<0.001 ) );
analysis.loadClosestNode([0.5 2.5],["ux" "uy"], [ -100 0 ] )
analysis.fixNodes( fixedEdgeSelector, ["ux" ] );
analysis.fixClosestNode([0 0],["uy"], 0);
analysis.plotCurrentLoad();
analysis.plotSupport();

analysis.printProblemInfo();

q = analysis.solve();

analysis.plotMaps(["uy" "ux" "exx" "sxx" "syy" "sxy" "sHM"],0.1);
fe.plotWired(mesh.nodes,q,0.1);
