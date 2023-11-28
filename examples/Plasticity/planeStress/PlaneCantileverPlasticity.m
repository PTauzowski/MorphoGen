clear;
close all;
res = 10;
l = 3;
tic
sfL4 = ShapeFunctionL16;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, 2*l, l, 2*res, res, sfL4.pattern);
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)))<0.001 );
loadedEdgeSelector = Selector( @(x)( abs(x(:,1) - 2*l)<0.001 ) );

fe=PlaneStressElastoPlasticElem( sfL4, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial( material );
fe.props.h=1;
fe.plot(mesh.nodes);
[I,J,V,Ksize] = fe.sparseMatrixAllocDataUniform( ["ux" "uy"] );
analysis = ElastoPlasticity( fe, mesh );
analysis.elementLoadLineIntegral( "global",loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + [0 150] ));
%supports = problem.fixNodes( fixedEdgeSelector, ["ux" ] ) | problem.fixClosestNode([0 0],["uy"], 0);
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysis.plotCurrentLoad();
analysis.plotSupport();
analysis.printProblemInfo();

analysis.solve();

analysis.plotMaps(["ux" "uy" "sxx" "syy" "sHM"],0.1);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);

