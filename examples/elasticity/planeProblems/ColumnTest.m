clear;
close all;
sf=ShapeFunctionL4;
l=1;
lw=1.4*l;
h=0.1;
b=0.1;
J=b*h^3/12;
nl=60;
E=205E9;
nu=0.3;
xp=[b/2 l];
Pcr=pi^2*E*J/lw^2
P=[0 -1];
nEigenForms=6;
%P=[0 1];

model = ColumnModel(sf,b,h,l,nl,E,nu,P,xp);
model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
figure;
stability=LinearStability(model.analysis.felems, model.mesh);
stability.Pfem=model.analysis.Pfem;
stability.supports=model.analysis.supports;
stability.solve(nEigenForms);
lambdas1=stability.lambdas;
lambdas1(1,1)
Pcr

for k=1:min(6,nEigenForms)
    figure;
    stability.setForm(k);
    model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);
    title(['Eigen form:' num2str(k)]);
end

analysis = SecondOrderElasticityWeighted( model.fe, model.mesh, false );
%analysis = LinearElasticityWeighted( model.fe, model.mesh, false );
analysis.Pnodal=model.analysis.Pnodal;
analysis.Pfem=model.analysis.Pfem;
analysis.supports=model.analysis.supports;

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


% tic
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.42, true );
% [objF, xopt]  = topOpt.solve();
% toc
% 
% figure;
% stability.solveWeighted(topOpt.x,5);
% lambdas2=stability.lambdas;
% stability.setForm(3);
% model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);