clear;
close all;


E=210E09;
nu=0.3;
R=0.25;
r=0.235;
alpha=22.5;
segmentLength=0.3;
res=15;
ShapeFn=ShapeFunctionL8;
frameElems=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
nArms=size(frameElems,1);

nSamples=5000;
samples=random("Uniform",0,360,nSamples,nArms);
samples(:,1)=0;

%[maxHM, endPoints, frameNodes] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, samples, ShapeFn);

load("ManipulatorOpti5000_2.mat");

tic;
[vN, vTz, vTy, vMs, vMz, vMy] = computeAllInternalForces(frameNodes,frameElems,E,nu);
toc

[vMinN1, iminN1]=min(vN(:,1));
[vMaxN1, imaxN1]=max(vN(:,1));
[vSortN1, iSortN1]=sort(vN(:,1));

[vMinN2, iminN2]=min(vN(:,2));
[vMaxN2, imaxN2]=max(vN(:,2));
[vSortN2, iSortN2]=sort(vN(:,2));


[vMinMs1, iminMs1]=min(vMs(:,1));
[vMaxMs1, imaxMs1]=max(vMs(:,1));
[vSortMs1, iSortMs1]=sort(vMs(:,1));

iminMs = iminMs1;

sampleMinMs = samples(iminMs1,:);
sampleMaxMs = samples(imaxMs1,:);

%save("ManipulatorOpti5000_2f.mat");

alpha=22.5;
r=0.24;
res=15;
Rfilter = 1.5*(R-r);
penal=3;
cutTreshold = 0.02;

[vMin, imin]=min(maxHM);
[vMax, imax]=max(maxHM);
[vSort, iSort]=sort(maxHM);

modelMin = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(imin,:), ShapeFn, true);
modelMax = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(iSort(nSamples),:), ShapeFn, true);

modelMinMs = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMinMs, ShapeFn, false);
modelMaxMs = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxMs, ShapeFn, false);


%manipulatorHMplot(modelMin,endPoints);
%manipulatorHMplot(modelMax,endPoints);

% nConfigs=5;
% plotExtremalConfigurations(E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSort,iSort);

plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(nSamples-0)),'MaxHM_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-0),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 2, \sigma_{HM}=" + num2str(vSort(nSamples-1)),'MaxHM_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-1),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 3, \sigma_{HM}=" + num2str(vSort(nSamples-2)),'MaxHM_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-2),:), ShapeFn, 0.0);

plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(1)),'MinHM_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(1),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 2, \sigma_{HM}=" + num2str(vSort(2)),'MinHM_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(2),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 3, \sigma_{HM}=" + num2str(vSort(3)),'MinHM_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(3),:), ShapeFn, 0.0);

plotArmConfigurationHMextended("maximal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(nSamples-0)),'MaxMs_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-0),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("maximal torsion moment configuration 2, M_{s}=" + num2str(vSortMs1(nSamples-1)),'MaxMs_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-1),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("maximal torsion moment configuration 3, M_{s}=" + num2str(vSortMs1(nSamples-2)),'MaxMs_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-2),:), ShapeFn, 0.0);

plotArmConfigurationHMextended("minimal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(1)),'MinMs_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(1),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("minimal torsion moment configuration 2, M_{s}=" + num2str(vSortMs1(2)),'MinMs_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(2),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("minimal torsion moment configuration 3, M_{s}=" + num2str(vSortMs1(3)),'MinMs_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(3),:), ShapeFn, 0.0);


plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(nSamples-0)),'BendingConfiguration.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-0),:), ShapeFn, 0.0);
plotArmConfigurationHMextended("maximal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(nSamples-0)),'TorsionConfiguration.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-0),:), ShapeFn, 0.0);

sampleMinA=[0 180 180 180 180 180 180]; %  0  164.4607  177.0448  215.1802  214.3430  186.1769  240.6389
sampleMaxA=[0 0 0 180 180 180 180];     %  0  343.7190   54.9144  125.4541  163.9858  169.3933  110.6999


% [maxHMa, endPointsa, frameNodes] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, [sampleMinA; sampleMaxA], ShapeFn);
% 
% plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMinA, ShapeFn);
% title(['Model for minimal [averaged] Huber-Mises for HMmax=' num2str(maxHMa(1))]);
% 
% plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxA, ShapeFn);
% title(['Model for maximal [averaged] Huber-Mises for HMmax=' num2str(maxHMa(2))]);

modelMinA = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMinA, ShapeFn, false);
modelMaxA = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxA, ShapeFn, false);


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

[xopt_bending, xopt_bending_buckling, lambda1, lambda2 ] = ArmTopOptBucklingFn('maxHM',R,r,segmentLength,alpha,Ty,Tz,N,My,Mz,Ms);

plotArmTopOptConfigProjections("BendingTopology","Bending topology",Rfilter, modelMaxMs.analysis, xopt_bending_buckling, cutTreshold, penal, false);

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

%load("ArmBendingTopOpt.mat");

% % xopt_full_nobuckling = plotArmTopOptConfig(Rfilter, modelMaxA.analysis, xopt, cutTreshold, penal, false);
% % title("Horizontal HM topology without buckling");
% % xopt_full_nobuckling = plotArmTopOptConfig(Rfilter, modelMaxA.analysis, xopt_buckling, cutTreshold, penal, false);
% % modelMaxA.analysis.mesh.exportMeshToFile(find(xopt_full_nobuckling>0.5),"OptTopologyArmMeshBending");
% % title("Horizontal HM topology with buckling");


loadFactor=0.5E8;

N  = vN(iminMs,1)*loadFactor;
Tz = vTz(iminMs,1)*loadFactor;
Ty = vTy(iminMs,1)*loadFactor;
Ms = vMs(iminMs,1)*loadFactor;
Mz = vMz(iminMs,1)*loadFactor;
My = vMy(iminMs,1)*loadFactor;

[xopt_torsion, xopt_torsion_buckling, lambda1, lambda2 ] = ArmTopOptBucklingFn('maxMs',R,r,segmentLength,alpha,Ty,Tz,N,My,Mz,Ms);

plotArmTopOptConfigProjections("TorsionTopology","Torsion topology",Rfilter, modelMaxMs.analysis, xopt_torsion_buckling, cutTreshold, penal, false);


% xopt_full_nobuckling = plotArmTopOptConfig(Rfilter, modelMaxMs.analysis, xopt, cutTreshold, penal, false);
% title("Max Ms topology without buckling");
% xopt_full_nobuckling = plotArmTopOptConfig(Rfilter, modelMaxMs.analysis, xopt_buckling, cutTreshold, penal, false);
% modelMaxA.analysis.mesh.exportMeshToFile(find(xopt_full_nobuckling>0.5),"OptTopologyArmMeshTorsion");
% title("Max Ms HM topology with buckling");




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


save("ManipulatorBucklingAndTorsion.mat");