clear;
close all;

% Cantilever topology optimization elastic task.

% Resolution of shortest (vertical) edge
res = 80;

% height of the cantilever
h = 1;

% Aspect ratio length/height
aspect=2;

% Filtering radius
Rfilter = 2*h/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

% Type of shape function to be used (here: four node Langrange)
sfL4 = ShapeFunctionL4;

% Creating FE mesh object
mesh = Mesh();

% Generating rectangular mesh ( aspect*h x h )
mesh.addRectMesh2D(0, 0, aspect*h, h, aspect*res, res, sfL4.pattern);

% Create plane stress finite element object
fe=PlaneStressElem( sfL4, mesh.elems );

% Create isotropic material object
material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);

% Assigning material to finite element
fe.setMaterial( material );

% Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
%analysis = LinearElasticityWeighted( fe, mesh, true );
analysis = SecondOrderElasticityWeighted(fe, mesh, 0, true);

% Creating node selector object to select fixed edge (left)
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.0005 ) );

% Fixing structure according to above defined node selector object
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );

% Creating load vector with one node loaded at the middle of right edge
analysis.loadClosestNode([aspect*h, h/2 ], ["ux" "uy"], [0 -0.0025] );

nEigenForms=10;
stability = LinearStability( analysis.felems, mesh);
stability.Pnodal = analysis.Pnodal;
stability.Pfem = analysis.Pfem;
stability.supports = analysis.supports;
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

%analysis = SecondOrderElasticityWeighted( fe, mesh, lambdas(nEigenForms), false );
analysis = SecondOrderElasticityWeighted( fe, mesh, 0.0, false );
analysis.Pnodal=stability.Pnodal;
analysis.Pfem=stability.Pfem;
analysis.supports=stability.supports;

figure;
tic
topOpt = StressIntensityTopologyOptimizationBuckling( Rfilter, analysis, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOpt.solve();
toc

% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
% [objF, xopt]  = topOpt.solve();
% toc
