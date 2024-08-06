clear;
close all;
sf=ShapeFunctionL9;
l=2;
h=5;
nl=200;
E=1;
nu=0.3;
P=[0 -1];
xp=[l/2 h];

model = ColumnModel(sf,l,h,nl,E,nu,P,xp);
model.analysis.solve(5);

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
figure;
model.analysis.setForm(5);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.2);

figure;
analysis = SecondOrderElasticityWeighted( model.fe, model.mesh, false );
analysis.Pfem=model.analysis.Pfem;

% Filtering radius
Rfilter = 2*h/nl;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;


figure;
tic
topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
[objF, xopt]  = topOpt.solve();
toc



tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, true );
[objF, xopt]  = topOpt.solve();
toc


