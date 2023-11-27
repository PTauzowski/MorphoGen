clear;
close all;
res = 20;
l = 3;
x0 = 5;
y0 = 5;
a=10;
hf = 0.5;
sf = ShapeFunctionL16;
mesh = Mesh();
mesh.addRectWithHoleMesh2D( 10, x0, y0, hf, res, sf.pattern );
%mesh.transformMeshDeg2D( [137 0], -90, [-137 0] );
fe=PlaneStressElem( sf, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(210000, 0.3);
fe.props.h=1;
fe.setMaterial( material );
fe.plot(mesh.nodes);

analysis = LinearElasticity( fe, mesh );
loadedEdge1 = Selector( @(x)( abs(x(:,1) -(x0 - a) ) < 0.001 ) );
loadedEdge2 = Selector( @(x)( abs(x(:,1) -(x0 + a) ) < 0.001 ) );
loadedEdge3 = Selector( @(x)( abs(x(:,2) -(y0 - a) ) < 0.001 ) );
loadedEdge4 = Selector( @(x)( abs(x(:,2) -(y0 + a) ) < 0.001 ) );
circleSelector = Selector( @(x)( abs(((x(:,1) - x0).^2 + (x(:,2) - y0).^2 ) - (a*hf)^2 )) < 0.01 );
analysis.elementLoadLineIntegral( "local", circleSelector,  ["ux" "uy"], @(x)( x*0 + [ 0 -100 ] ));
analysis.fixClosestNode( [ x0-a y0-a], ["ux" "uy"], [0 0] );
analysis.fixClosestNode([ x0+a y0-a],"uy", 0 );
analysis.plotCurrentLoad();
%problem.plotSupport();

analysis.printProblemInfo();

analysis.solve();
figure
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);
analysis.plotMaps(["uy" "ux" "sxx" "syy" "s1" "s2" "sxy" "sHM"],0.1);
