clear;
close all;
sf=ShapeFunctionL27;
l=1;
b=0.1;
h=0.1;
lw=2*l;
nl=80;
E=205E9;
nu=0.3;
rho=7850; 
xp=[b/2 h/2 l];
Pcr=pi^2*E*b*h^3/12/lw^2;
P=[0 0 -1];
nEigenForms=20;
%P=[0 1];

model = ColumnModel3D( sf,b,h,l,nl,E,nu,rho,P,xp);
model.solveWeighted();

% model.plotModel();
% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "szz" "sHM"],0.1);
% model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
stability = LinearStability( model.analysis.felems, model.mesh);
stability.Pnodal = model.analysis.Pnodal;
stability.Pfem = model.analysis.Pfem;
stability.supports = model.analysis.supports;
stability.solve( nEigenForms);
lambdas1 = stability.lambdas;
disp(['Eigenvalue = ', num2str(lambdas1(1,1))]),
disp(['P_critical = ', num2str(Pcr)])

for k=1:min( 6, nEigenForms)
    figure;
    stability.setForm(k);
    model.fe.plotWired( model.mesh.nodes,stability.qnodal, 0.2);
    title(['Eigenvector: ' num2str(k)]);
end

% vibrations=LinearNaturalVibration(model.analysis.felems, model.mesh);
% vibrations.Pnodal=model.analysis.Pnodal;
% vibrations.Pfem=model.analysis.Pfem;
% vibrations.supports=model.analysis.supports;
% vibrations.solve(nEigenForms);
% omegas=sqrt(vibrations.lambdas);
% omegas(1,1)
% 
% 
% for k=1:min(6,nEigenForms)
%     figure;
%     vibrations.setForm(k);
%     model.fe.plotSolidDeformed(model.mesh.nodes,vibrations.qnodal,0.2);
%     title(['Eigen form:' num2str(k)]);
% end

vibrations=LinearNaturalVibration(model.analysis.felems, model.mesh);
vibrations.Pnodal=model.analysis.Pnodal;
vibrations.Pfem=model.analysis.Pfem;
vibrations.supports=model.analysis.supports;
vibrations.solve(nEigenForms);
omegas=sqrt(vibrations.lambdas);
3.5160/l^2*sqrt(E*J)
omegas(1,1)


for k=1:min(6,nEigenForms)
    figure;
    vibrations.setForm(k);
    model.fe.plotSolidDeformed(model.mesh.nodes,vibrations.qnodal,0.2);
    title(['Eigen form:' num2str(k)]);
end

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