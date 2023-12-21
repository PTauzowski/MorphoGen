clear;
close all;

% Rectangular beam 
% resolution on height
res = 40;

% beam length
l = 1;

% filter diameter
Rfilter = 4*l/res;

%penalty factor
penalty = 3;


% stress treshold ratio for element's removal
cutTreshold=0.005;

% * FEM task definition *
tic
aspect=3;
sfL4 = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, aspect*l, l, aspect*res, res, sfL4.pattern);
fe=PlaneStressElem( sfL4, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);
fe.setMaterial( material );
fe.props.h=1;

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1))<0.001 ) );
loadedEdgeSelector = Selector( @(x)( (abs( x(:,1) - aspect*l) < 1.0E-4) & ( abs(x(:,2)-0)<l/15 )  ) );
analysis.loadClosestNode([0, l ], ["ux" "uy"], [0 -1] );
analysis.fixNodes( fixedEdgeSelector, [ "ux" ] );
analysis.fixClosestNode([aspect*l 0], ["uy"], 0) ;

analysis.printProblemInfo();
fe.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();

tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penalty, 0.4, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penalty, 0.4, true);
[objF, xopt]  = topOpt.solve();
toc




