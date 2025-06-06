clear;
clear classes
clear all
rehash
close all;

% Cantilever topology optimization elastic task.

% Resolution of shortest (vertical) edge
res = 40;

% height of the cantilever
h = 1;

% Aspect ratio length/height
aspect=6;

% length of the cantilever
l=aspect*h;

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
mesh.addRectMesh2D(0, 0, l, h, aspect*res, res, sfL4.pattern);

% Create plane stress finite element object
fe=PlaneStressElem( sfL4, mesh.elems );

% Create isotropic material object
material = PlaneStressMaterial('mat1');
material.setElasticIzo(205E9, 0.3);
material.setMassIzoMatrix(7850);
% Assigning material to finite element
fe.setMaterial( material );

% Ke = fe.computeStifnessMatrix(mesh.nodes,1,1);
% Me = fe.computeMassMatrix(mesh.nodes,1,1);
% 
%  [q, l] = eigs( Ke(:,:,1),Me(:,:,1), 8, 'smallestabs');

% Creating linear elastic finite element analysis object with weighted matrix feature, weighted by element density.
analysisLinear = LinearElasticityWeighted( fe, mesh, false );

% Creating node selector object to select fixed edge (left)
loadEdgeSelector = Selector( @(x)( abs(x(:,2) -h ) < 0.0005 ) );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.0005 ) );
fixedEdgeSelector2 = Selector( @(x)( abs(x(:,1)-l) < 0.0005 ) );

% Fixing structure according to above defined node selector object
analysisLinear.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
analysisLinear.fixNodes( fixedEdgeSelector2,  "ux" );

% Creating load vector with one node loaded at the middle of right edge
P=-2.0E9; %100;
%P=-1.5E8; %150;

hp=h;
analysisLinear.loadClosestNode([ l, hp ], ["ux" "uy"], [0 P] );
%analysisLinear.elementLoadLineIntegral( "global",loadEdgeSelector, ["ux" "uy"], @(x)( x*0 + [0 P/l] ));


const_rows=3;
ncel=round(const_rows*res*l);
const_elems=[1:ncel size(mesh.elems,1):-1:size(mesh.elems,1)-ncel ];


nEigenForms=10;
stability = LinearStability( analysisLinear.felems, mesh);
stability.Pnodal = analysisLinear.Pnodal;
stability.Pfem = analysisLinear.Pfem;
stability.supports = analysisLinear.supports;
% stability.solve( nEigenForms);
% lambdas = diag(stability.lambdas);
% for k=1:min(10,nEigenForms)
%     subplot(5, 2, k);
%     stability.setForm(k);
%     fe.plotWired(mesh.nodes,stability.qnodal,0.2);
%     axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
%     lambda_str = sprintf('%.4g', lambdas(k));
%     title(['Form:' num2str(k), ' \lambda=' lambda_str]);
% end

nEigenForms=10;
vibrations = LinearNaturalVibration( analysisLinear.felems, mesh );
vibrations.Pnodal = analysisLinear.Pnodal;
vibrations.Pfem = analysisLinear.Pfem;
vibrations.supports = analysisLinear.supports;
vibrations.solve( nEigenForms);
omegas = diag(vibrations.omegas)/2/pi;
% for k=1:min(10,nEigenForms)
%     subplot(5, 2, k);
%     %figure
%     vibrations.setForm(1,k);
%     fe.plotWithSettings(mesh.nodes,"deformed",vibrations.qnodal,0.1);
%     %axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
%     omega_str = sprintf('%.4g', omegas(k));
%     title(['Form:' num2str(k), ' \omega=' omega_str]);
% end

Pin1 = omegas(1)^2*vibrations.fromFEMVector(vibrations.qforms(:,1));
Pin2 = omegas(2)^2*vibrations.fromFEMVector(vibrations.qforms(:,2));
Pin3 = omegas(3)^2*vibrations.fromFEMVector(vibrations.qforms(:,3));
Pin5 = omegas(5)^2*vibrations.fromFEMVector(vibrations.qforms(:,5));
Pin7 = omegas(7)^2*vibrations.fromFEMVector(vibrations.qforms(:,7));
Pin9 = omegas(9)^2*vibrations.fromFEMVector(vibrations.qforms(:,9));
Pin10 = omegas(10)^2*vibrations.fromFEMVector(vibrations.qforms(:,10));


harmonicVivrations = ElasticHarmonicVibrations(fe, mesh, 1, false, true);
harmonicVivrations.Pnodal=analysisLinear.Pnodal;
harmonicVivrations.Pfem=analysisLinear.Pfem;
harmonicVivrations.supports=analysisLinear.supports;

alphas = [0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6];

% conds1 = harmonicVivrations.tabMatrixCondition(alphas,1);
% conds2 = harmonicVivrations.tabMatrixCondition(alphas,2);
% conds3 = harmonicVivrations.tabMatrixCondition(alphas,3);

% load('matrixConditions.mat');
% 
% figure, hold on
% p1=plot(alphas,conds1,'b','LineWidth', 3);
% title('Matrix K-\alpha\omega_1^2*M condition ');
% xlabel('\alpha');
% ylabel('rcond(K-\alpha\omega_1^2*M )');
% %xlim([37 57]);
% set(gca, 'FontSize', 16)
% 
% figure, hold on
% p1=plot(alphas,conds2,'b','LineWidth', 3);
% title('Matrix K-\alpha\omega_2^2*M condition ');
% xlabel('\alpha');
% ylabel('rcond(K-\alpha\omega_2^2*M )');
% %xlim([37 57]);
% set(gca, 'FontSize', 16)
% 
% figure, hold on
% p1=plot(alphas,conds3,'b','LineWidth', 3);
% title('Matrix K-\alpha\omega_3^2*M condition ');
% xlabel('\alpha');
% ylabel('rcond(K-\alpha\omega_3^2*M )');
% %xlim([37 57]);
% set(gca, 'FontSize', 16)

% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc
% 
% analysisLinear.Pnodal = Pin1;
% analysisLinear.Pfem=analysisLinear.toFEMVector(Pin1);
% 
% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc
% 
% analysisLinear.Pnodal = Pin2;
% analysisLinear.Pfem=analysisLinear.toFEMVector(Pin2);
% 
% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc
% 
% analysisLinear.Pnodal = Pin3;
% analysisLinear.Pfem=analysisLinear.toFEMVector(Pin3);
% 
% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc
% 
% analysisLinear.Pnodal = Pin5;
% analysisLinear.Pfem=analysisLinear.toFEMVector(Pin5);
% 
% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc
% 
% analysisLinear.Pnodal = Pin7;
% analysisLinear.Pfem=analysisLinear.toFEMVector(Pin7);
% 
% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc

% analysisLinear.Pnodal = Pin9;
% analysisLinear.Pfem=analysisLinear.toFEMVector(Pin9);
% 
% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc
% 
% analysisLinear.Pnodal = Pin10;
% analysisLinear.Pfem=analysisLinear.toFEMVector(Pin10);
% 
% figure;
% tic
% topOptLinear = StressIntensityTopologyOptimizationVol( Rfilter, analysisLinear, cutTreshold, penal, 0.5, true );
% topOptLinear.setConstElems(const_elems);
% [objF, xopt]  = topOptLinear.solve();
% toc

% 
% figure;
% tic
% topOptHarmonicVivrations = StressIntensityTopologyOptimizationDynamicBuckling( Rfilter, harmonicVivrations, cutTreshold, penal, 0.45, false );
% topOptHarmonicVivrations.setConstElems(const_elems);
% toc

basename='Const_load';
nforms=3;
for k=1:nforms
    harmonicVivrations.count=0;
    harmonicVivrations.mode=k;
    topOptHarmonicVivrations = StressIntensityTopologyOptimizationDynamicBuckling( Rfilter, harmonicVivrations, cutTreshold, penal, 0.45, false );
    topOptHarmonicVivrations.setConstElems(const_elems);
    [objF, xopt]  = topOptHarmonicVivrations.solve();
    topOptHarmonicVivrations.plot_frequencies(basename, nforms)
    topOptHarmonicVivrations.plot_forms(basename, nforms, topOptHarmonicVivrations.iteration-1)
    topOptHarmonicVivrations.plot_correlation_curve( basename, nforms)
    topOptHarmonicVivrations.plot_uncorrelated_frames( basename)
    topOptHarmonicVivrations.resetAnalysis();
    harmonicVivrations.modes=[];
end

harmonicVivrations.isLoadConst=false;

basename='Var_load';
for k=1:nforms
    harmonicVivrations.count=0;
    harmonicVivrations.mode=k;
    topOptHarmonicVivrations = StressIntensityTopologyOptimizationDynamicBuckling( Rfilter, harmonicVivrations, cutTreshold, penal, 0.45, false );
    topOptHarmonicVivrations.setConstElems(const_elems);
    [objF, xopt]  = topOptHarmonicVivrations.solve();
    topOptHarmonicVivrations.plot_frequencies(basename, nforms)
    topOptHarmonicVivrations.plot_forms(basename, nforms, topOptHarmonicVivrations.iteration-1)
    topOptHarmonicVivrations.plot_correlation_curve( basename, nforms)
    topOptHarmonicVivrations.plot_uncorrelated_frames( basename)
    topOptHarmonicVivrations.resetAnalysis();
    harmonicVivrations.modes=[];
end

% [K, M] = harmonicVivrations.computeMatrices(xopt);
% 
% save KM_matrices.mat K M

% figure;
% tic
% topOptSecondOrder = StressIntensityTopologyOptimizationDynamicBuckling( Rfilter, analysisSecondOrder, cutTreshold, penal, 0.30, true );
% topOptSecondOrder.setConstElems(const_elems);
% [objF, xopt]  = topOptSecondOrder.solve();
% toc

% figure;
% tic
% topOptBuckling = StressIntensityTopologyOptimizationBuckling( Rfilter, analysisWithBuckling, cutTreshold, penal, 0.38, true );
% [objF, xopt]  = topOptBuckling.solve();
% toc

% figure;
% U = normalize(analysisSecondOrder.modes2, 1); 
% C = abs(U' * U);
% imagesc(C);
% colormap(flipud(hot));
% colorbar;
% axis equal tight;
% xlabel('Mode number');
% ylabel('Mode number');
% title('Correlation matrix of eigenmode 1');
% 
% figure;
% U = normalize(analysisSecondOrder.modes2, 1); 
% C = abs(U' * U);
% imagesc(C);
% colormap(flipud(hot));
% colorbar;
% axis equal tight;
% xlabel('Mode number');
% ylabel('Mode number');
% title('Correlation matrix of eigenmode 2');
% 
% figure;
% U = normalize(analysisSecondOrder.modes3, 1); 
% C = abs(U' * U);
% imagesc(C);
% colormap(flipud(hot));
% colorbar;
% axis equal tight;
% xlabel('Mode number');
% ylabel('Mode number');
% title('Correlation matrix of eigenmode 3');



%save('Cantilever2DBucklingDown80.mat');

