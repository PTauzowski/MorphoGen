clear;
close all;
res = 15;
l = 3;
tic
ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addRectMesh3D( 0, 0, 0, 2*l, l, l, 2*res, res, res, ShapeFn.localNodes);
fe = SolidElasticElem( ShapeFn, mesh.elems );
mX = max(mesh.nodes(:,1));

fe.plotSolid(mesh.nodes);
fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material)

analysis = LinearElasticity( fe, mesh );

fixedFaceSelector = Selector( @(x)( x(:,1) ) );
loadedFaceSelector = Selector( @(x)( x(:,1) - 2*l ) );
analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy"], @(x)( x(:,1:2)*0 + [0 -100] ));
analysis.fixNodes( fixedFaceSelector, ["ux" "uy" "uz"] );

analysis.printProblemInfo();
analysis.plotCurrentLoad();
analysis.plotSupport();

tic
q = analysis.solve();
toc

analysis.plotMaps([ "sxx" "szz" "sHM"],0.1);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);


