clear;
close all;
res = 20;
l = 2;
h = 1;
E=210000;
nu=0.3;
sy=190;
tic
sfL4 = ShapeFunctionL9;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, l, h, round(l/h*res), res, sfL4.pattern);
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)))<0.001 );
loadedEdgeSelector = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );

fe=PlaneStressElastoPlasticElem( sfL4, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElastoPlasticIzo(E, nu, sy);
fe.setMaterial( material );
fe.plot(mesh.nodes);
%[I,J,V,Ksize] = fe.sparseMatrixllocDataUniform( ["ux" "uy"] );
analysis = ElastoPlasticity( fe, mesh );
analysis.elementLoadLineIntegral( "global",loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + [0 -45] ));
%supports = problem.fixNodes( fixedEdgeSelector, ["ux" ] ) | problem.fixClosestNode([0 0],["uy"], 0);
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysis.plotCurrentLoad();
analysis.plotSupport();
analysis.printProblemInfo();

analysis.solve();

analysis.plotMaps(["ux" "uy" "sxx" "syy" "sHM" "pz"],0.1);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);

