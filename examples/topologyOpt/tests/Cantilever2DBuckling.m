clear;
close all;

% Cantilever topology optimization elastic task.

% Resolution of shortest (vertical) edge
res = 50;

% height of the cantilever
h = 1;

% Aspect ratio length/height
aspect=2;

% Filtering radius
Rfilter = 4*h/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

% Type of shape function to be used (here: four node Langrange)
sfL4 = ShapeFunctionL4;

% Creating FE mesh object
mesh = Mesh();

% Generating rectangular mesh ( aspect*h x h )
elems = mesh.addRectMesh2D(0, 0, aspect*h, h, aspect*res, res, sfL4.pattern);

% Create plane stress finite element object
fe=PlaneStressElem( sfL4, elems );

% Create isotropic material object
material = PlaneStressMaterial('mat1');
material.setElasticIzo(205E9, 0.3);

% Assigning material to finite element
fe.setMaterial( material );

% Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
analysisLinear      = LinearElasticityWeighted( fe, mesh, false );
analysisSecondOrder = SecondOrderElasticityWeighted(fe, mesh, 0, false);

% Creating node selector object to select fixed edge (left)
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.0005 ) );

% Fixing structure according to above defined node selector object
analysisLinear.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysisSecondOrder.fixNodes( fixedEdgeSelector, ["ux" "uy"] );

% Creating load vector with one node loaded at the middle of right edge
P=-2.0E8; %100;
%P=-1.5E8; %150;

hp=0;
analysisLinear.loadClosestNode([aspect*h, hp ], ["ux" "uy"], [0 P] );
analysisSecondOrder.loadClosestNode([aspect*h, hp ], ["ux" "uy"], [0 P] );

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


analysisWithBuckling = SecondOrderElasticityWeighted( fe, mesh, 0.95, false );
analysisWithBuckling.Pnodal=stability.Pnodal;
analysisWithBuckling.Pfem=stability.Pfem;
analysisWithBuckling.supports=stability.supports;

figure;
tic
topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOptLinear.solve();
toc

figure;
tic
topOptSecondOrder = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.38, true );
[objF, xopt]  = topOptSecondOrder.solve();
toc

figure;
tic
topOptBuckling = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisWithBuckling, cutTreshold, penal, 0.38, true );
[objF, xopt]  = topOptBuckling.solve();
toc

%load('Cantilever2DBucklingDown80.mat');

figure, hold on
p1=plot(topOptSecondOrder.plVol,topOptSecondOrder.plLambda,'b','LineWidth', 3);
set(gca, 'XDir', 'reverse');
title('Critical force coefficient evolution without buckling');
xlabel('Volume fracion [%]');
ylabel('Critical force coefficient [%]');
xlim([37 57]);
set(gca, 'FontSize', 24)

figure, hold on
p2=plot(topOptBuckling.plVol,topOptBuckling.plLambda,'r','LineWidth', 3);
set(gca, 'XDir', 'reverse');
title('Critical force coefficient evolution with buckling');
xlabel('Volume fracion [%]');
ylabel('Critical force coefficient [%]');
xlim([37 57]);
set(gca, 'FontSize', 24)

figure, hold on
p1=plot(topOptSecondOrder.plVol,topOptSecondOrder.plLambda,'b','LineWidth', 3);
p2=plot(topOptBuckling.plVol,topOptBuckling.plLambda,'r','LineWidth', 3);
legend([p1, p2], {'Without buckling', 'With buckling'});
set(gca, 'XDir', 'reverse');
title('Comparison of critical force evolution coefficient');
xlabel('Volume fracion [%]');
ylabel('Critical force coefficient [%]');
xlim([37 57]);
set(gca, 'FontSize', 24)

% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
% [objF, xopt]  = topOpt.solve();
% toc

save('Cantilever2DBucklingDown80.mat');

