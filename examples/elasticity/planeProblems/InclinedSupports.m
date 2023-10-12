clear;
close all;
res = 20;
l = 3;
segmentHeigh=l;
segmentLength=2*l;
tic;
endAngle=45;
segmentShortening=segmentHeigh/2*tan(endAngle*pi/180);
sfL4 = ShapeFunctionL4;
mesh = Mesh();
%mesh.addRectMesh2D(0, 0, 2*l, l, 2*res, res, sfL4.pattern);
mesh.addRectMesh2D(0, -segmentHeigh/2, segmentLength, segmentHeigh, 2*res, res, sfL4.pattern);
mesh.transformNodesXY( @(x)([ x(:,1)+2*segmentShortening.*x(:,2)./segmentHeigh.*(x(:,1))./segmentLength  x(:,2) ]) )
fe=PlaneStressElem( sfL4, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial( material );
fe.props.h=1;
fe.plotSolid(mesh.nodes);

analysis = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( @(x)( x(:,1) ) );
inclinedSupportEdgeSelector = Selector( @(x)( x(:,1)-(segmentLength+2*segmentShortening.*x(:,2)./segmentHeigh ) ) );
loadedEdgeSelector = Selector( @(x)( x(:,2) - segmentHeigh/2 ) );
circleSelector = Selector( @(x)( ((x(:,1) - 1.5).^2 + (x(:,2) - 1.5).^2 )-1.5 ), 0.2 );
%fedges = fe.findEdges( problem.findNodes( loadedEdgeSelector ));
%problem.plotSelectedNodeNumbers( loadedEdgeSelector );
analysis.elementLoadLineIntegral("global",loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + [0 -150] ));
%supports = problem.fixNodes( fixedEdgeSelector, ["ux" ] ) | problem.fixClosestNode([0 0],["uy"], 0);
%supports = problem.fixNodes( fixedEdgeSelector, ["ux" "uy"] ) | problem.fixClosestNode([segmentLength 0],["uy"], 0) ;
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysis.fixNodes( inclinedSupportEdgeSelector, ["uy"] );
analysis.setNodalSupportRotations( inclinedSupportEdgeSelector, 90-endAngle)
analysis.plotCurrentLoad();
analysis.plotSupport();
analysis.printProblemInfo();

q = analysis.solve();

analysis.plotMaps(["ux" "uy" "sxx" "syy" "sHM"],0.0);
fe.plotWired(mesh.nodes,analysis.qnodal,0.1);


