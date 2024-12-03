clear;
close all;

E=210E09;
nu=0.3;
R=0.25;
r=0.20;
segmentLength=0.3;
res=10;
alpha=30;
betas=[ 23 18 70 35 ];

% Filtering radius
Rfilter = 1.5*(R-r);
penal=3;
cutTreshold = 0.05;

ShapeFn = ShapeFunctionL8;
model = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, betas, ShapeFn);

analysis = LinearElasticityWeighted( model.fe, model.mesh, false );
%problem = LinearElasticity( fe, mesh );
fixedEdgeSelector = Selector( model.fixedSurfaceNodes );
loadedFaceSelector = Selector( model.loadSurfaceNodes );
%constElemsSelector =  Selector( @(x)( (x(:,3) < 0.05 * Length ) ) & (x(:,3) > 0.96 * Length ) );

analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 -100] ));
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);
const_elems = model.const_elems;

%model.fe.plot(model.mesh.nodes);
%problem.plotNodes();
%analysis.plotCurrentLoad();
%analysis.plotSupport();

view(45, 45);

analysis.printProblemInfo();

tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, true);
[objF, xopt]  = topOpt.solve();
toc




