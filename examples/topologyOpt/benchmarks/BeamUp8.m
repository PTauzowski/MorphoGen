clear;
close all;
res = 80;
l = 1;
aspect=8;
lp=1;

Rfilter = 2*l/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

sfL4 = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, aspect*l, l, aspect*res, res, sfL4.pattern);
fe=PlaneStressElem( sfL4, mesh.elems );
material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);
fe.setMaterial( material );
fe.props.h=1;
P=0.25;

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.001 ) );
loadedEdgeSelector = Selector( @(x)( (abs( x(:,1) - aspect*l) < 1.0E-4) & ( abs(x(:,2)-0)<l/15 ) ) );
%analysis.loadClosestNode([0 l ], ["ux" "uy"], [P -P] );
analysis.loadClosestNode([lp l ], ["ux" "uy"], [0 -P] );
analysis.loadClosestNode([2*lp, l ], ["ux" "uy"], [0 -P] );
analysis.loadClosestNode([3*lp, l ], ["ux" "uy"], [0 -P] );
analysis.loadClosestNode([4*lp, l ], ["ux" "uy"], [0 -P] );
analysis.loadClosestNode([5*lp, l ], ["ux" "uy"], [0 -P] );
analysis.loadClosestNode([6*lp, l ], ["ux" "uy"], [0 -P] );
analysis.loadClosestNode([7*lp, l ], ["ux" "uy"], [0 -P] );
%analysis.loadClosestNode([8*lp, l ], ["ux" "uy"], [-P -P] );
analysis.fixClosestNode([0 0], ["ux" "uy"], [0 0]);
analysis.fixClosestNode([aspect*l 0], ["uy"], 0) ;

analysis.printProblemInfo();
fe.plot(mesh.nodes);
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

