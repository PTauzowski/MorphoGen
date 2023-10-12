clear;
close all;
res = 10;
l = 3;
tic
sfT3 = ShapeFunctionT3;
mesh = Mesh();
P = [0 0; 1 0; 1 0.7; 2 2; 2 3; 1 3; 1 2; 0 1];
C = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1];
%xv = [1 4 4 1 1 NaN 2 2 3 3 2];
%yv = [1 1 4 4 1 NaN 2 3 3 2 2];
mesh.addDelaunayMesh2D( P, C, 500 );
fe=PlaneStressElem( sfT3, mesh.elems );
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
%problem.elementLoadLineIntegral( "global", loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + [ -100 0 ] ));
analysis.loadClosestNode([0.5 2.5],["ux" "uy"], [ -100 0 ] )
analysis.fixNodes( fixedEdgeSelector, ["ux" ] );
analysis.fixClosestNode([0 0],["uy"], 0);
analysis.plotCurrentLoad();
analysis.plotSupport();

analysis.printProblemInfo();

q = analysis.solve();

analysis.plotMaps(["uy" "ux" "exx" "sxx" "syy" "sxy" "sHM"],0.1);
fe.plotWired(mesh.nodes,q,0.1);
