clear;
close all;

res = 50;
l = 1;
tic
aspect=3;

Rfilter = 2*l/res;

% stress treshold ratio for element's removal
cutTreshold=0.005;


%penalty factor
penal = 3;


sfL4 = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, aspect*l, l, aspect*res, res, sfL4.pattern);
fe=PlaneStressElem( sfL4, mesh.elems );

material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.0);
fe.setMaterial( material );

fe.props.h=1;
analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( x(:,1) ) );
loadedEdgeSelector = Selector( @(x)( not( (abs( x(:,1) - aspect*l) < 1.0E-4) & ( abs(x(:,2)-0)<l/15 ) ) ) );
analysis.loadClosestNode([aspect*l/2, 0 ], ["ux" "uy"], [0 -1] );
analysis.fixClosestNode([0 0], ["ux" "uy"], [0 0]);
analysis.fixClosestNode([aspect*l 0], ["uy"], [0]) ;
analysis.printProblemInfo();
fe.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();

tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, true);
[objF, xopt]  = topOpt.solve();
toc





