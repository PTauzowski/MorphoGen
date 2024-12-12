clear;
close all;


E=210E09;
nu=0.3;
R=0.25;
r=0.24;
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

r=0.24;
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


frameElems=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
mesh=Mesh();
frameElem=Frame3D(frameElems,E,0.02,0.8*E,0.0004,0.0004,0.003);
mesh.nodes=modelMin.frameNodes;
[Fel, Feg] = computeInternalForces(frameElem,mesh);

lb=[0 0 0 0 0 0];
ub=[360 360 360 360 360 360];
x0=[180 180 180 180 180 180];

xopt = optimizeArmHM(E,nu,segmentLength,R,r,res, alpha, lb, ub, x0, ShapeFn);

disp("Global internal forces");
printInternalForcesTable(Feg);

disp("Local internal forces");
printInternalForcesTable(Fel);

barNumber=2;
loadFactor=1.0E8;

N  = Fel(1,barNumber)*loadFactor;
Ty = Fel(2,barNumber)*loadFactor;
Tz = Fel(3,barNumber)*loadFactor;
Ms = Fel(4,barNumber)*loadFactor;
My = Fel(5,barNumber)*loadFactor;
Mz = Fel(6,barNumber)*loadFactor;

%ArmTopOptBucklingFn('maximal1',R,r,2*segmentLength,alpha,Ty,Tz,N,My,Mz,Ms);
 
% Filtering radius
Rfilter = 1.5*(R-r);
penal=3;
cutTreshold = 0.02;

% tic
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, modelMax.analysis, cutTreshold, penal, 0.4, false );
% topOpt.setConstElems(modelMax.const_elems);
% [objF, xopt]  = topOpt.solve();
% 
% 
% topOpt.plotCurrentFrame();
% toc
% 
% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, modelMax.analysis, penal, 0.3, false);
% [objF, xopt]  = topOpt.solve();
% toc
% 
% 


