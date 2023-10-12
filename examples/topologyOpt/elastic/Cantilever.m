clear;
close all;
res = 100;
l = 1;

aspect=2;

Rfilter = 3*l/res;

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

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( x(:,1) ) );
loadedEdgeSelector = Selector( @(x)( not( (abs( x(:,1) - aspect*l) < 1.0E-4) & ( abs(x(:,2)-l)<l/15 ) ) ) );
analysis.loadClosestNode([aspect*l, l/2 ], ["ux" "uy"], [0 -1] );
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );

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



