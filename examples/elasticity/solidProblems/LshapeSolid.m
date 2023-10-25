clear;
close all;
res = 10;
l = 3;
tic
ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addLshape3D( 2*l, 0.8*l, res, ShapeFn.localNodes );
fixedFaceSelector = Selector( @(x)( abs(x(:,3) - 2*l)<0.001 ) );
loadedFaceSelector = Selector( @(x)( abs(x(:,1) - 2*l)<0.001 ) );

fe = SolidElasticElem( ShapeFn, mesh.elems );
mX = max(mesh.nodes(:,1));

fe.plotSolid(mesh.nodes);
fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material)

analysis = LinearElasticity( fe, mesh );

analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy"], @(x)( x(:,1:2)*0 + [0 -100] ));
analysis.fixNodes( fixedFaceSelector, ["ux" "uy" "uz"] );

analysis.printProblemInfo();
analysis.plotCurrentLoad();
analysis.plotSupport();

tic
q = analysis.solve();
toc

analysis.plotMaps([ "sxx" "szz" "sHM"],0.3);
fe.plotWired(mesh.nodes,analysis.qnodal,0.3);



