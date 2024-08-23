clear;
close all;
sf=ShapeFunctionL4;
l=1;
lw=2*l;
h=0.1;
b=0.1;
J=b*h^3/12;
nl=300;
E=205E9;
nu=0.3;
xp=[b/2 l];
Pcr=pi^2*E*J/lw^2
P=[0 1];
nEigenForms=6;

model = ColumnModel(sf,b,h,l,nl,E,nu,P,xp);
% model.solveWeighted();
% model.plotModel();
% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
% model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

stability=LinearStability(model.analysis.felems, model.mesh);
stability.Pnodal=model.analysis.Pnodal;
stability.Pfem=model.analysis.Pfem;
stability.supports=model.analysis.supports;
stability.solve(nEigenForms);
lambdas=diag(stability.lambdas)
disp(['Eigenvalue = ', num2str(lambdas(1))]),
disp(['P_critical = ', num2str(Pcr)])

P_critical_str = sprintf('%.4g', Pcr);
sgtitle(['P_{cr} = ', P_critical_str]);
for k=1:min(6,nEigenForms)
    subplot(3, 2, k);
    stability.setForm(k);
    model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);
    lambda_str = sprintf('%.4g', lambdas(k));
    title(['Form:' num2str(k), ' \lambda=' lambda_str]);
end

analysis = SecondOrderElasticityWeighted( model.fe, model.mesh, lambdas(nEigenForms), false );
%analysis = LinearElasticityWeighted( model.fe, model.mesh, false );
analysis.Pnodal=model.analysis.Pnodal;
analysis.Pfem=model.analysis.Pfem;
analysis.supports=model.analysis.supports;

% model.solveWeighted();
% 
% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
% model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

% Filtering radius
Rfilter = 2*h/nl;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

Pb=[0 1];
modelb = ColumnModel(sf,b,h,l,nl,E,nu,Pb,xp);

% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.2, true);
% [objF, xopt]  = topOpt.solve();
% toc


tic
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.3, true );
[objF, xopt]  = topOpt.solve();
toc

figure;
stability.solveWeighted(topOpt.x,6);
lambdas2=diag(stability.lambdas);

disp(['Eigenvalue = ', num2str(lambdas2(1))]),
disp(['P_critical = ', num2str(Pcr)])

P_critical_str = sprintf('%.4g', Pcr);
sgtitle(['P_{cr} = ', P_critical_str]);
for k=1:min(6,nEigenForms)
    subplot(3, 2, k);
    stability.setForm(k);
    model.fe.plotWired(model.mesh.nodes,stability.qnodal,0.2);
    lambda_str = sprintf('%.4g', lambdas(k));
    title(['Form:' num2str(k), ' \lambda=' lambda_str]);
end