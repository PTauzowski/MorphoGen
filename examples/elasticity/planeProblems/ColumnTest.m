clear;
close all;
sf=ShapeFunctionL4;
l=1;
h=5;
nl=200;
E=2.1E9;
nu=0.3;
P=[0 -10000];
xp=[l/2 h];

model = ColumnModel(sf,l,h,nl,E,nu,P,xp);
model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
figure;
stability=LinearStability(model.analysis.felems, model.mesh);
stability.Pfem=model.analysis.Pfem;
stability.solve(5);
lambdas1=stability.lambdas;
stability.setForm(3);
model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);

figure;
analysis = SecondOrderElasticityWeighted( model.fe, model.mesh, false );
analysis = LinearElasticityWeighted( model.fe, model.mesh, false );
analysis.Pnodal=model.analysis.Pnodal;

% Filtering radius
Rfilter = 2*h/nl;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;


% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
% [objF, xopt]  = topOpt.solve();
% toc


tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
stability.solveWeighted(topOpt.x,5);
lambdas2=stability.lambdas;
stability.setForm(3);
model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);