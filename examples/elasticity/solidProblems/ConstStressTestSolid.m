clear;
close all;
res = 10;
l = 3;
tic
ShapeFn = ShapeFunctionL8;
%ShapeFn = ShapeFunctionT4;
mesh = Mesh();
mesh.addRectMesh3D( 0, 0, 0, 2.5*l, l, l, 2.5*res, res, res, ShapeFn.localNodes);
%mesh.addRectMeshTetrahedral3D( '6T', [0 0 0], [2.5*l l l], [2.5*res, res, res] )
fe = SolidElasticElem( ShapeFn, mesh.elems );
mX = max(mesh.nodes(:,1));

%fe.plotSolid(mesh.nodes);

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material)

analysis = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1) ) < 0.0001 ) );
loadedFaceSelector = Selector( @(x)( abs(x(:,1) - 2.5*l) < 0.0001 ) );
analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [-2 0 0] ));
analysis.fixNodes( fixedEdgeSelector, ["ux"] );
analysis.fixClosestNode([0 0 0],["uy" "uz"], 0);

analysis.printProblemInfo();
%problem.plotCurrentLoad();
%problem.plotSupport();


q  = analysis.solve();


analysis.plotMaps([ "sxx" "szz" "sHM"],0.0);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);

