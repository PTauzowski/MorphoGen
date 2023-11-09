clear;
close all;

res = 20;
l = 1;

Rfilter = 1.2*l/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

sfL4 = ShapeFunctionL4;
mesh = Mesh();
mesh.addLshape(2*l, 0.8*l, 2*res, sfL4.pattern);
fe=PlaneStressElem( sfL4, mesh.elems );

material = PlaneStressMaterial('mat1');
fe.props.h=1;
material.setElasticIzo(1, 0.3);
fe.setMaterial( material );

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( abs(x(:,2) - 2*l) < 0.001 ) );
loadedEdgeSelector = Selector( @(x)( (abs( x(:,1) - 2*l) < 1.0E-4) & ( abs(x(:,2)-0.4*l)<l/30 ) ) );

analysis.loadClosestNode([ 2*l, 0.4*l ], ["ux" "uy"], [0 -1.0] );
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

