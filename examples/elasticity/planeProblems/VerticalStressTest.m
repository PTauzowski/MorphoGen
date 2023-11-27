clear;
%close all;
res = 2;
l = 3;
tic
sfL4 = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, l, 2*l, res, 2*res, sfL4.pattern);
fe=PlaneStressElem( sfL4, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.setMaterial( material );
fe.props.h=1;
fe.plot(mesh.nodes);
problem = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( @(x)( abs(x(:,2))<0.001 ) );
loadedEdgeSelector = Selector( @(x)( abs(x(:,2) - 2*l) <0.001) );
%problem.plotSelectedNodeNumbers( loadedEdgeSelector );
problem.elementLoadLineIntegral( "global", loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + [ 0 -100 ] ));
problem.fixNodes( fixedEdgeSelector, ["uy" ] );
problem.fixClosestNode([0 0],["ux"], 0);
problem.plotCurrentLoad();
problem.plotSupport();

problem.printProblemInfo();

q = problem.solve();

problem.plotMaps(["ux" "uy" "sxx" "syy" "sHM"],0.2);
fe.plotWired(mesh.nodes,problem.qnodal,0.1);
