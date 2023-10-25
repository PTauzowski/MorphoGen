clear;
close all;
l=3;
c=0.4;
res = 24;
sf = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, l, c*l, round(1/c*res), res, sf.pattern);
mesh.addRectMesh2D(0, c*l, c*l, (1-c)*l, res, 1/(1-c)*res, sf.pattern);
fe=PlaneStressElem( sf, mesh.elems );
fe.props.h=1;
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.setMaterial( material );

fe.plotSolid(mesh.nodes);
analysis = LinearElasticity( fe, mesh );
fixedEdgeSelectorX = Selector( @(x)( abs(x(:,1))<0.001 ) );
fixedEdgeSelectorY = Selector( @(x)( abs(x(:,2)) )<0.001 );
loadedEdgeSelectorX = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );
loadedEdgeSelectorY = Selector( @(x)( abs(x(:,2) - l)<0.001 ) );
%problem.plotSelectedNodeNumbers( loadedEdgeSelector );
analysis.elementLoadLineIntegral( "global", loadedEdgeSelectorX, ["ux" "uy"], @(x)( x*0 + [ 100 0 ] ));
analysis.elementLoadLineIntegral( "global", loadedEdgeSelectorY, ["ux" "uy"], @(x)( x*0 + [ 0 100 ] ));
analysis.fixNodes( fixedEdgeSelectorX, "ux");
analysis.fixNodes( fixedEdgeSelectorY, "uy" );
analysis.plotCurrentLoad();
analysis.plotSupport();

analysis.printProblemInfo();
analysis.solve();

analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);
