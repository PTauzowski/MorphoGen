clear;
close all;
sf=ShapeFunctionL8;
l=10;
b=1;
h=1;
lw=2*l;
nl=100;
E=1;
nu=0.3;
xp=[b/2 h/2 l];
Pcr=pi^2*E*b*h^3/12/lw^2
P=[0 0 -0.8*Pcr/h/b];
%P=[0 1];

model = ColumnModel3D(sf,b,l,nl,E,nu,P,xp);
model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "szz" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
figure;
stability=LinearStability(model.analysis.felems, model.mesh);
stability.Pfem=model.analysis.Pfem;
stability.supports=model.analysis.supports;
stability.solve(5);
lambdas1=stability.lambdas
stability.setForm(1);
model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);

analysis = SecondOrderElasticityWeighted( model.fe, model.mesh, false );
%analysis = LinearElasticityWeighted( model.fe, model.mesh, false );
analysis.Pnodal=model.analysis.Pnodal;
analysis.Pfem=model.analysis.Pfem;
analysis.supports=model.analysis.supports;

% Filtering radius
Rfilter = 2*l/nl;

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
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.42, true );
[objF, xopt]  = topOpt.solve();
toc
% 
% figure;
% stability.solveWeighted(topOpt.x,5);
% lambdas2=stability.lambdas;
% stability.setForm(3);
% model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);