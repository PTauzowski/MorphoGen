clear;
close all;
res = 30;
l = 3;

% Filtering radius
Rfilter = 2*l/res;
penal=3;
cutTreshold = 0.005;

ShapeFn = ShapeFunctionL8;
mesh = Mesh();
mesh.addRectMesh3D( 0, 0, 0, 2*l, l/2, l, 2*res, res/2, res, ShapeFn.localNodes);
fe = SolidElasticElem( ShapeFn, mesh.elems );

fe.props.h=1;
material = SolidMaterial('mat1');
material.setElasticIzo(1, 0.3);
material.setElasticIzoGrad();
fe.setMaterial(material)

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( (x(:,1)==0) & (x(:,3)==0) ) );
lengthSymmetrySelector = Selector( @(x)( abs(x(:,1) - 2*l) < 0.001 ) );
depthSymmetrySelector = Selector( @(x)( abs(x(:,2) - l/2) < 0.001 ) );
analysis.loadClosestNode( [2*l 0 0], "uz", -1 );
analysis.loadClosestNode( [2*l 0 l], "uy", 0.5 );
analysis.loadClosestNode( [l 0 0], "uz", -1 );
analysis.loadClosestNode( [l 0 l], "uy", 0.5 );
analysis.fixNodes( fixedEdgeSelector, ["uz"] );
analysis.fixNodes( lengthSymmetrySelector, ["ux"] );
analysis.fixNodes( depthSymmetrySelector, ["uy"] );

fe.plotSolid(mesh.nodes);
%problem.plotNodes();
analysis.plotCurrentLoad();
analysis.plotSupport();

view(45, 45);

analysis.printProblemInfo();

tic
%topOpt = StressIntensityTopologyOptimizationVol( Rfilter, problem, cutTreshold, penal, 0.4, true );
%[objF, xopt]  = topOpt.solve();
toc

figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
[objF, xopt]  = topOpt.solve();
toc

