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

[vMinTz1, iminTz1]=min(vTz(:,1));
[vMaxTz1, imaxTz1]=max(vTz(:,1));
[vSortTz1, iSortTz1]=sort(vTz(:,1));

[vMinTy1, iminTy1]=min(vTy(:,2));
[vMaxTy1, imaxTy1]=max(vTy(:,2));
[vSortTy1, iSortTy1]=sort(vTy(:,2));

iminMs = iminMs1;

sampleMinMs = samples(iminMs1,:);
sampleMaxMs = samples(imaxMs1,:);

sampleMinTz = samples(iminTz1,:);
sampleMaxTz = samples(imaxTz1,:);

sampleMinTy = samples(iminTy1,:);
sampleMaxTy = samples(imaxTy1,:);

%save("ManipulatorOpti5000_2f.mat");

alpha=22.5;
r=0.22;
res=15;
Rfilter = 1.5*(R-r);
penal=3;
cutTreshold = 0.02;

[vMin, imin]=min(maxHM);
[vMax, imax]=max(maxHM);
[vSort, iSort]=sort(maxHM);

modelMin = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(imin,:), ShapeFn, true);
modelMax = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(iSort(nSamples),:), ShapeFn, true);

modelMinTz = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMinTz, ShapeFn, false);
modelMaxTz = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxTz, ShapeFn, false);

modelMinTy = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMinTy, ShapeFn, false);
modelMaxTy = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxTy, ShapeFn, false);


%manipulatorHMplot(modelMin,endPoints);
%manipulatorHMplot(modelMax,endPoints);

nConfigs=5;

%plotExtremalConfigurations(E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSortMs1,iSortMs1);

% plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(nSamples-0)),'MaxHM_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-0),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 2, \sigma_{HM}=" + num2str(vSort(nSamples-1)),'MaxHM_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 3, \sigma_{HM}=" + num2str(vSort(nSamples-2)),'MaxHM_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-2),:), ShapeFn, 0.0);
% 
% plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(1)),'MinHM_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 2, \sigma_{HM}=" + num2str(vSort(2)),'MinHM_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(2),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal Huber-Mises stress configuration 3, \sigma_{HM}=" + num2str(vSort(3)),'MinHM_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(3),:), ShapeFn, 0.0);
% 
% plotArmConfigurationHMextended("maximal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(nSamples-0)),'MaxMs_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-0),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal torsion moment configuration 2, M_{s}=" + num2str(vSortMs1(nSamples-1)),'MaxMs_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal torsion moment configuration 3, M_{s}=" + num2str(vSortMs1(nSamples-2)),'MaxMs_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-2),:), ShapeFn, 0.0);
% 
% plotArmConfigurationHMextended("minimal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(1)),'MinMs_configuration1.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(1),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal torsion moment configuration 2, M_{s}=" + num2str(vSortMs1(2)),'MinMs_configuration2.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(2),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("minimal torsion moment configuration 3, M_{s}=" + num2str(vSortMs1(3)),'MinMs_configuration3.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(3),:), ShapeFn, 0.0);
% 
% 
% plotArmConfigurationHMextended("maximal Huber-Mises stress configuration 1, \sigma_{HM}=" + num2str(vSort(nSamples-0)),'BendingConfiguration.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSort(nSamples-0),:), ShapeFn, 0.0);
% plotArmConfigurationHMextended("maximal torsion moment configuration 1, M_{s}=" + num2str(vSortMs1(nSamples-0)),'TorsionConfiguration.pdf',E,nu,segmentLength,R,r,res, modelMax.halfSegmentNelems, alpha, samples(iSortMs1(nSamples-0),:), ShapeFn, 0.0);

sampleMinN=[0 180 180 180 180 180 180]; %  0  164.4607  177.0448  215.1802  214.3430  186.1769  240.6389
sampleMaxMz=[0 0 0 180 180 180 180];     %  0  343.7190   54.9144  125.4541  163.9858  169.3933  110.6999
sampleMaxTy=[ 0  0  180 0 180  180  180]; %  0    0.0253  246.2603  287.9922  340.8601  106.8641  112.8238
sampleMaxMs=[ 0  45 45 45  270  180 180]; %  0 0   33.5628   56.5082   36.2258  278.0525  195.3981  122.2441

[maxHMa, endPointsa, frameNodes] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, [ sampleMaxMz; sampleMaxTy; sampleMaxMs], ShapeFn);
% 
plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxMz, ShapeFn);
title(['Model for minimal [averaged] Huber-Mises for HMmax=' num2str(maxHMa(1))]);

plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxTy, ShapeFn);
title(['Model for maximal shear force Ty max=' num2str(maxHMa(2))]);

plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxMs, ShapeFn);
title(['Model for maximal torsion moment Ms max=' num2str(maxHMa(3))]);

loadFactor=0.5E8;
[xopt_bending, xopt_bending_buckling, lambda1, lambda2 ]  = configurationTopology(E,nu,R,r,segmentLength,ShapeFn,alpha,sampleMaxMz,frameElems,loadFactor);
[xopt_shear, xopt_shear_buckling, lambda1, lambda2 ]  = configurationTopology(E,nu,R,r,segmentLength,ShapeFn,alpha,sampleMaxTy,frameElems,loadFactor);
[xopt_torsion, xopt_torsion_buckling, lambda1, lambda2 ]  = configurationTopology(E,nu,R,r,segmentLength,ShapeFn,alpha,sampleMaxMs,frameElems,loadFactor);

modelMz = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxMz, ShapeFn, false);
modelTy = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxTy, ShapeFn, false);
modelMs = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sampleMaxMs, ShapeFn, false);

plotArmTopOptConfigProjections("BendingTopology","Bending topology",Rfilter, modelMz.analysis, xopt_bending_buckling, cutTreshold, penal, false);
plotArmTopOptConfigProjections("ShearTopology","Shear topology",Rfilter, modelTy.analysis, xopt_shear_bucling, cutTreshold, penal, false);
plotArmTopOptConfigProjections("TorsionTopology","Bending topology",Rfilter, modelMs.analysis, xopt_torsion_buckling, cutTreshold, penal, false);






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


%save("ManipulatorBucklingAndTorsion.mat");