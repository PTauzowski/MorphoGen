clear;
close all;


E=210E09;
nu=0.3;
R=0.25;
r=0.23;
segmentLength=0.3;
res=15;
alpha=30;
ShapeFn=ShapeFunctionL8;
nArms=7;


nSamples=5000;
samples=random("Uniform",0,360,nSamples,nArms);
samples(:,1)=0;

%[maxHM, endPoints] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, samples, ShapeFn);

load("ManipulatorOpti5000.mat");

r=0.23;
res=10;

[vMin, imin]=min(maxHM);
[vMax, imax]=max(maxHM);
[vSort, iSort]=sort(maxHM);

modelMin = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(imin,:), ShapeFn);
modelMax = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(iSort(nSamples-2),:), ShapeFn);

%manipulatorHMplot(modelMin,endPoints);
%manipulatorHMplot(modelMax,endPoints);

nConfigs=2;
plotExtremalConfigurations(E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSort,iSort);

sampleMinA=[0 180 180 180 180 180 180]; %  0  164.4607  177.0448  215.1802  214.3430  186.1769  240.6389
sampleMaxA=[0 0 0 180 180 180 180];     %  0  343.7190   54.9144  125.4541  163.9858  169.3933  110.6999


[maxHMa, endPointsa] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, [sampleMinA; sampleMaxA], ShapeFn);

plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMinA, ShapeFn);
title(['Model for minimal [averaged] Huber-Mises for HMmax=' num2str(maxHMa(1))]);

plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxA, ShapeFn);
title(['Model for maximal [averaged] Huber-Mises for HMmax=' num2str(maxHMa(2))]);

modelMinA = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMinA, ShapeFn);
modelMaxA = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxA, ShapeFn);

frameElems=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
mesh=Mesh();
frameElem=Frame3D(frameElems,E,0.02,0.8*E,0.0004,0.0004,0.003);
mesh.nodes=modelMaxA.frameNodes;
[Fel, Feg] = computeInternalForces(frameElem,mesh);

lb=[0 0 0 0 0 0];
ub=[360 360 360 360 360 360];
x0=[180 180 180 180 180 180];

%xopt = optimizeArmHM(E,nu,segmentLength,R,r,res, alpha, lb, ub, x0, ShapeFn);

disp("Global internal forces");
printInternalForcesTable(Feg);

disp("Local internal forces");
printInternalForcesTable(Fel);

barNumber=2;
loadFactor=0.5E8;

N  = Fel(7,barNumber)*loadFactor;
Tz = Fel(8,barNumber)*loadFactor;
Ty = Fel(9,barNumber)*loadFactor;
Ms = Fel(10,barNumber)*loadFactor;
Mz = Fel(11,barNumber)*loadFactor;
My = Fel(12,barNumber)*loadFactor;

% [xopt, xopt_buckling, lambda1, lambda2 ] = ArmTopOptBucklingFn('maximal1',R,r,segmentLength,alpha,Ty,Tz,N,My,Mz,Ms);
% 
% % Filtering radius
% Rfilter = 1.5*(R-r);
% penal=3;
% cutTreshold = 0.02;
% 
% plotArmTopOptConfig(Rfilter, modelMin.analysis, xopt, cutTreshold, penal, false);
% title("Minimal HM topology without buckling");
% plotArmTopOptConfig(Rfilter, modelMin.analysis, xopt_buckling, cutTreshold, penal, false);
% title("Minimal HM topology with buckling");
% 
% plotArmTopOptConfig(Rfilter, modelMax.analysis, xopt, cutTreshold, penal, false);
% title("Maximal HM topology without buckling");
% plotArmTopOptConfig(Rfilter, modelMax.analysis, xopt_buckling, cutTreshold, penal, false);
% title("Maximal HM topology with buckling");
% 
% plotArmTopOptConfig(Rfilter, modelMinA.analysis, xopt, cutTreshold, penal, false);
% title("Vertical HM topology without buckling");
% plotArmTopOptConfig(Rfilter, modelMinA.analysis, xopt_buckling, cutTreshold, penal, false);
% title("Vertical HM topology with buckling");
% 
% plotArmTopOptConfig(Rfilter, modelMaxA.analysis, xopt, cutTreshold, penal, false);
% title("Horizontal HM topology without buckling");
% plotArmTopOptConfig(Rfilter, modelMaxA.analysis, xopt_buckling, cutTreshold, penal, false);
% title("Horizontal HM topology with buckling");



%[objF, xopt]  = topOpt.solve();


% toc
% 
% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, modelMax.analysis, penal, 0.3, false);
% [objF, xopt]  = topOpt.solve();
% toc
% 
% 


