clear;
close all;
res = 60;
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
P=0.00020;

% Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
analysisLinear = LinearElasticityWeighted( fe, mesh, false );


fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.001 ) );
loadedEdgeSelector = Selector( @(x)( (abs( x(:,1) - aspect*l) < 1.0E-4) & ( abs(x(:,2)-0)<l/15 ) ) );
%analysis.loadClosestNode([0 l ], ["ux" "uy"], [P -P] );
analysisLinear.loadClosestNode([lp l ], ["ux" "uy"], [0 -P] );
analysisLinear.loadClosestNode([2*lp, l ], ["ux" "uy"], [0 -P] );
analysisLinear.loadClosestNode([3*lp, l ], ["ux" "uy"], [0 -P] );
analysisLinear.loadClosestNode([4*lp, l ], ["ux" "uy"], [0 -P] );
analysisLinear.loadClosestNode([5*lp, l ], ["ux" "uy"], [0 -P] );
analysisLinear.loadClosestNode([6*lp, l ], ["ux" "uy"], [0 -P] );
analysisLinear.loadClosestNode([7*lp, l ], ["ux" "uy"], [0 -P] );
%analysis.loadClosestNode([8*lp, l ], ["ux" "uy"], [-P -P] );
analysisLinear.fixClosestNode([0 0], ["ux" "uy"], [0 0]);
analysisLinear.fixClosestNode([aspect*l 0], ["uy"], 0) ;

analysisLinear.printProblemInfo();
fe.plot(mesh.nodes);
analysisLinear.plotCurrentLoad();
analysisLinear.plotSupport();


nEigenForms=10;
stability = LinearStability( analysisLinear.felems, mesh);
stability.Pnodal = analysisLinear.Pnodal;
stability.Pfem = analysisLinear.Pfem;
stability.supports = analysisLinear.supports;
stability.solve( nEigenForms);
lambdas = diag(stability.lambdas)
for k=1:min(10,nEigenForms)
    subplot(5, 2, k);
    stability.setForm(k);
    fe.plotWired(mesh.nodes,stability.qnodal,0.2);
    axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
    lambda_str = sprintf('%.4g', lambdas(k));
    title(['Form:' num2str(k), ' \lambda=' lambda_str]);
end

analysisWithBuckling = SecondOrderElasticityWeighted( fe, mesh, 0.90, false );
analysisWithBuckling.Pnodal=stability.Pnodal;
analysisWithBuckling.Pfem=stability.Pfem;
analysisWithBuckling.supports=stability.supports;


analysisSecondOrder = SecondOrderElasticityWeighted(fe, mesh, 0, false);
analysisSecondOrder.Pnodal=stability.Pnodal;
analysisSecondOrder.Pfem=stability.Pfem;
analysisSecondOrder.supports=stability.supports;



% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.4, true );
% [objF, xopt]  = topOptLinear.solve();
% toc

% figure;
% tic
% topOptSecondOrder = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.4, true );
% [objF, xopt]  = topOptSecondOrder.solve();
% toc

figure;
tic
topOptBuckling = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisWithBuckling, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOptBuckling.solve();
toc




% tic
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.4, true );
% [objF, xopt]  = topOpt.solve();
% toc
% 
% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysisLinear, penal, 0.4, true);
% [objF, xopt]  = topOpt.solve();
% toc
% 
