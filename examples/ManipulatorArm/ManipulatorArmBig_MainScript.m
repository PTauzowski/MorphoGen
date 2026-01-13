clear;
close all;


% PLA (printed)
E = 2.0e9; % Pa
nu = 0.35;
rho = 1250; % kg/m³

ShapeFn=ShapeFunctionL8;
frameElems=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
nArms=size(frameElems,1);
% 
% nSamples=100000;
% samples=random("Uniform",0,360,nSamples,nArms);
% samples(:,1)=0;

nDiv  = 8;                           % 8 bins per angle
%vals  = linspace(0, 360 - 360/nDiv, nDiv);  % [0,45,90,...,315]
vals  = linspace(-180, 180 - 360/nDiv, nDiv);  % [0,45,90,...,315]

% We keep the first column fixed at 0, vary the remaining 6
nVar = nArms - 1;                    % 6
G = cell(1, nVar);
[G{:}] = ndgrid(vals);

% Assemble samples: rows = 8^6, cols = 7 (first col fixed to 0)
samples = zeros(nDiv^nVar, nArms);
%samples(:,1) = 0;
for k = 1:nVar
    samples(:, k+1) = G{k}(:);
end

nSamples = size(samples,1);

% frameNodes = modelRef.computeFrameNodesBatch( segmentLength, alpha, samples );

R=0.14; % m, outer radius
r=0.105; % m, inner radius
alpha=22.5; % deg, segment connection inclination
segmentLength=0.15; % m, segment length
res=15;
res_thickness=3;

Pz = 100; % N -  Vertical force at the end of manipulator;

%[~, ~, frameNodes] = computeArmSamples(E, nu, segmentLength, alpha, samples );

%save("ManipulatorFrameConfigurations200K.mat","frameNodes");
load("ManipulatorFrameConfigurations200K.mat","frameNodes");

max(max(max(frameNodes)))

% figure;
% hold on, axis on; 
% daspect([1 1 1]);
% xlabel("x");
% ylabel("y");
% zlabel("z");
% scatter3(endPoints(:,1),endPoints(:,2),endPoints(:,3),'Marker','.');
% scatter3(endPoints(:,1),endPoints(:,2),endPoints(:,3)*0-3,'Marker','.','MarkerEdgeColor','r');

modelRef  = ManipulatorModel3D(E, nu, segmentLength, R, r, res, res_thickness, alpha, samples(1,:), ShapeFn, true, Pz);
modelRef.fe.plot(modelRef.mesh.nodes);
daspect([1 1 1]);

%[vN, vTy, vTz, vMs, vMy, vMz, all_forces] = computeAllInternalForces( frameNodes, frameElems, E, nu, R, r, samples);

%save("ManipulatorFrameForces200K.mat","vN", "vTy", "vTz", "vMs", "vMy", "vMz", "all_forces");
load("ManipulatorFrameForces200K.mat","vN", "vTy", "vTz", "vMs", "vMy", "vMz", "all_forces");

%[vMaxN1, imaxN1]=max(abs(vN(:,1)));
%[vMaxN2, imaxN2]=max(abs(vN(:,2)));

maxLoad=zeros(6,2);

[max_N,  max_conf_N,  max_seg_N  ] = findMaxConfiguration(squeeze(vN(:,1,:)));
[max_Ms, max_conf_Ms, max_seg_Ms ] = findMaxConfiguration(squeeze(vMs(:,1,:)));
[max_Ty, max_conf_Ty, max_seg_Ty ] = findMaxConfiguration(squeeze(vTy(:,1,:)));
[max_Tz, max_conf_Tz, max_seg_Tz ] = findMaxConfiguration(squeeze(vTz(:,1,:)));
[max_My, max_conf_My, max_seg_My ] = findMaxConfiguration(squeeze(vMy(:,1,:)));
[max_Mz, max_conf_Mz, max_seg_Mz ] = findMaxConfiguration(squeeze(vMz(:,1,:)));

gen_forces_maxN = squeeze(all_forces(max_conf_N,1,max_seg_N,1:6))';
gen_forces_maxTy = squeeze(all_forces(max_conf_Ty,1,max_seg_Ty,1:6))';
gen_forces_maxTz = squeeze(all_forces(max_conf_Tz,1,max_seg_Tz,1:6))';
gen_forces_maxMs = squeeze(all_forces(max_conf_Ms,1,max_seg_Ms,1:6))';
gen_forces_maxMy = squeeze(all_forces(max_conf_My,1,max_seg_My,1:6))';
gen_forces_maxMz = squeeze(all_forces(max_conf_Mz,1,max_seg_Mz,1:6))';

config_forces1 = squeeze(all_forces(max_conf_N,1,:,:))';
config_forces2 = squeeze(all_forces(max_conf_Ty,1,:,:))';
config_forces3 = squeeze(all_forces(max_conf_Tz,1,:,:))';
config_forces4 = squeeze(all_forces(max_conf_Ms,1,:,:))';
config_forces5 = squeeze(all_forces(max_conf_My,1,:,:))';
config_forces6 = squeeze(all_forces(max_conf_Mz,1,:,:))';

save("config_forces.mat",'config_forces1','config_forces2','config_forces3', 'config_forces4','config_forces5','config_forces6');

max_segments = [ max_seg_N max_seg_Ty max_seg_Tz max_seg_Ms max_seg_My max_seg_Mz ];

% [vMaxN1, imaxN1]=max(abs(vN(:,1,:)));
% [vSortN1, iSortN1]=sort(abs(vN(:,1,:)));
% 
% [vMaxMs1, imaxMs1]=max(abs(vMs(:,1,:)));
% [vSortMs1, iSortMs1]=sort(abs(vMs(:,1,:)));
% 
% [vMaxTz1, imaxTz1]=max(abs(vTz(:,1,:)));
% [vSortTz1, iSortTz1]=sort(abs(vTz(:,1,:)));
% 
% [vMaxTy1, imaxTy1]=max(abs(vTy(:,1,:)));
% [vSortTy1, iSortTy1]=sort(vTy(:,1,:));
% 
% [vMaxMz1, imaxMz1]=max(abs(vMz(:,1,:)));
% [vSortMz1, iSortMz1]=sort(abs(vMz(:,1,:)));
% 
% [vMaxMy1, imaxMy1]=max(abs(vMy(:,1,:)));
% [vSortMy1, iSortMy1]=sort(vMy(:,1,:));

sampleMaxN  = samples(max_conf_N,:);
sampleMaxTy = samples(max_conf_Ty,:);
sampleMaxTz = samples(max_conf_Tz,:);
sampleMaxMs = samples(max_conf_Ms,:);
sampleMaxMy = samples(max_conf_My,:);
sampleMaxMz = samples(max_conf_Mz,:);

% [vMin, imin]=min(maxHM);
% [vMax, imax]=max(maxHM);
% [vSort, iSort]=sort(maxHM);

modelMaxN  = ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, sampleMaxN,  ShapeFn, true, Pz);
modelMaxTy = ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, sampleMaxTy, ShapeFn, true, Pz);
modelMaxTz = ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, sampleMaxTz, ShapeFn, true, Pz);
modelMaxMs = ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, sampleMaxMs, ShapeFn, true, Pz);
modelMaxMy = ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, sampleMaxMy, ShapeFn, true, Pz);
modelMaxMz = ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, sampleMaxMz, ShapeFn, true, Pz);

%  modelMaxN.plotConfigurations('max_Q1', 'Configuration for extremal Q_1', gen_forces_maxN,sampleMaxN,max_seg_N);
%  modelMaxN.plotMesh();
% % 
% modelMaxN.compute(1);
% modelMaxN.plotHM_map(1)
% % 
% modelMaxTy.plotConfigurations('max_Q2','Configuration for extremal Q_2',gen_forces_maxTy,sampleMaxTy,max_seg_Ty);
% modelMaxTz.plotConfigurations('max_Q3','Configuration for extremal Q_3',gen_forces_maxTz,sampleMaxTz,max_seg_Tz);
% modelMaxMs.plotConfigurations('max_Q4','Configuration for extremal Q_4',gen_forces_maxMs,sampleMaxMs,max_seg_Ms);
% modelMaxMy.plotConfigurations('max_Q5','Configuration for extremal Q_5',gen_forces_maxMy,sampleMaxMy,max_seg_My);
% modelMaxMz.plotConfigurations('max_Q6','Configuration for extremal Q_6',gen_forces_maxMz,sampleMaxMz,max_seg_Mz);

modelMaxTy.saveModelMatrices("max_Q2");
% 
% fig = figure;
% hold on, axis on; 
% daspect([1 1 1]);
% xlabel("x");
% ylabel("y");
% zlabel("z");
% light('Position', [-1 -2 5], 'Style', 'local');
% light('Position', [1 1 5], 'Style', 'infinite');
% gca.FontSize = 18;

% modelMaxTy.plotMesh();
% modelMaxTz.plotMesh();
% modelMaxMs.plotMesh();
% modelMaxMy.plotMesh();
% modelMaxMz.plotMesh();
% 
% scatter3(endPointsNew(:,1),endPointsNew(:,2),endPointsNew(:,3),'Marker','.');
% scatter3(endPointsNew(:,1),endPointsNew(:,2),endPointsNew(:,3)*0-3,'Marker','.','MarkerEdgeColor',[0.9 0.9 0.9]);
% scatter3(endPointsNew(:,1),endPointsNew(:,2)*0+3,endPointsNew(:,3),'Marker','.','MarkerEdgeColor',[0.9 0.9 0.9]);
% scatter3(endPointsNew(:,1)*0-3,endPointsNew(:,2),endPointsNew(:,3),'Marker','.','MarkerEdgeColor',[0.9 0.9 0.9]);

% for k=1:20
%     model  = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(k+100,:),  ShapeFn, true);
%     model.plotMesh();
% end

% view(3);
% exportgraphics(fig, "Multiconfig.png", 'Resolution', 1200); 
% savefig(fig, "Multiconfig.fig");

%[objF, xopt] = computeFullArmTopOpt(E,nu,segmentLength,R,r,res, alpha, sampleMaxTz, ShapeFn, Rfilter, cutTreshold, penal);

%manipulatorHMplot(modelMax,endPoints);

% nConfigs=1;
% plotExtremalConfigurations("Extremal N", E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSortN1,iSortN1);
% plotExtremalConfigurations("Extremal T_z", E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSortTz1,iSortTz1);
% plotExtremalConfigurations("Extremal T_y", E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSortTy1,iSortTy1);
% plotExtremalConfigurations("Extremal M_y", E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSortMy1,iSortMy1);
% plotExtremalConfigurations("Extremal M_z", E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSortMz1,iSortMz1);
% plotExtremalConfigurations("Extremal M_s", E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSortMs1,iSortMs1);

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

sampleMinN_smooth=[0 180 180 180 180 180 180]; %  0  164.4607  177.0448  215.1802  214.3430  186.1769  240.6389
sampleMaxMz_smooth=[0 0 0 180 180 180 180];     %  0  343.7190   54.9144  125.4541  163.9858  169.3933  110.6999
sampleMaxTy_smooth=[ 0  0  180 0 180  180  180]; %  0    0.0253  246.2603  287.9922  340.8601  106.8641  112.8238
sampleMaxMs_smooth=[ 0  45 45 45  270  180 180]; %  0 0   33.5628   56.5082   36.2258  278.0525  195.3981  122.2441

% sampleMinN=[ 0  164.4607  177.0448  215.1802  214.3430  186.1769  240.6389 ];
% sampleMaxMz=[  0  343.7190   54.9144  125.4541  163.9858  169.3933  110.6999 ];
% sampleMaxTy=[ 0    0.0253  246.2603  287.9922  340.8601  106.8641  112.8238 ];
% sampleMaxMs=[  0   33.5628   56.5082   36.2258  278.0525  195.3981  122.2441 ];

modeSamples = [sampleMaxN; sampleMaxTy; sampleMaxTz; sampleMaxMz; sampleMaxMy; sampleMaxMs; -sampleMaxN; -sampleMaxTy; -sampleMaxTz; -sampleMaxMz; -sampleMaxMy; -sampleMaxMs];
real_max_segments = [max_segments max_segments];

%[maxHMa, endPointsa, frameNodes] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, modeSamples, ShapeFn);

% plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxMz, ShapeFn);
% title(['Model for minimal [averaged] Huber-Mises for HMmax=' num2str(vMaxMs1)]);
% 
% plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxTy, ShapeFn);
% title(['Model for maximal shear force Ty max=' num2str(vMaxTz1)]);
% 
% plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sampleMaxMs, ShapeFn);
% title(['Model for maximal torsion moment Ms max=' num2str(vMaxMs1)]);

% loadFactor=0.4E8;
% [xopt_bending, xopt_bending_buckling, bending_linear_lambda1, bending_buckling_lambda2 ]  = configurationTopology(E,nu,R,r,segmentLength,ShapeFn,alpha,sampleMaxMz,frameElems,loadFactor);
% [xopt_shear, xopt_shear_buckling, shear_linear_lambda1, shear_buckling_lambda2 ]  = configurationTopology(E,nu,R,r,segmentLength,ShapeFn,alpha,sampleMaxTy,frameElems,loadFactor);
% [xopt_torsion, xopt_torsion_buckling, torsion_linear_lambda1, torsion_buckling_lambda2 ]  = configurationTopology(E,nu,R,r,segmentLength,ShapeFn,alpha,sampleMaxMs,frameElems,loadFactor);

% fprintf('\n');
% disp(['Bending lambda linear = ' num2str(bending_linear_lambda1) ' Bending lambda buckling = ' num2str(bending_buckling_lambda2)]);
% disp(['Shear lambda linear   = ' num2str(shear_linear_lambda1) ' Shear lambda buckling   = ' num2str(shear_buckling_lambda2)]);
% disp(['Torsion lambda linear = ' num2str(torsion_linear_lambda1) ' Bending lambda buckling = ' num2str(torsion_buckling_lambda2)]);

E=2.0E9;
nu=0.4;

% plotArmTopOptConfigProjections("BendingTopology","Bending topology",Rfilter, modelMz.analysis, modelMz.halfSegmentNelems, 2, xopt_bending_buckling, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("ShearTopology","Shear topology",Rfilter, modelTy.analysis, modelTy.halfSegmentNelems, 2, xopt_shear_buckling, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("TorsionTopology","Torsion topology",Rfilter, modelMs.analysis, modelMs.halfSegmentNelems, 2, xopt_torsion_buckling, cutTreshold, penal, false);
%plotArmTopOptConfigProjections("FinalTopology","Final topology",Rfilter, modelMz.analysis, modelMz.halfSegmentNelems, 2, xopt_final, cutTreshold, penal, false);

%xOnes = xopt_bending_buckling;
xOnes(:)=1;
segmentNo=2;
dispFactor=0.0;
loadFactor=1.0;

modelMaxN.analysis.felems{1}.edge_color='k';
modelMaxN.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxN.analysis.felems{1}.face_alpha=1;
modelMaxTy.analysis.felems{1}.edge_color='k';
modelMaxTy.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxTy.analysis.felems{1}.face_alpha=1;
modelMaxTz.analysis.felems{1}.edge_color='k';
modelMaxTz.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxTz.analysis.felems{1}.face_alpha=1;
modelMaxMy.analysis.felems{1}.edge_color='k';
modelMaxMy.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxMy.analysis.felems{1}.face_alpha=1;
modelMaxMz.analysis.felems{1}.edge_color='k';
modelMaxMz.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxMz.analysis.felems{1}.face_alpha=1;
modelMaxMs.analysis.felems{1}.edge_color='k';
modelMaxMs.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxMs.analysis.felems{1}.face_alpha=1;

% plotArmTopOptConfigProjections("NormalInitialConfiguration","Normal initial configuration N",Rfilter, modelMaxN.analysis, modelMaxN.halfSegmentNelems, max_segments(1), 1, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("ShearInitialConfiguration","Shear initial configuration T_y",Rfilter, modelMaxTy.analysis, modelMaxTy.halfSegmentNelems, max_segments(2), 1, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("ShearInitialConfiguration","Shear initial configuration T_z",Rfilter, modelMaxTz.analysis, modelMaxTz.halfSegmentNelems, max_segments(3), 1, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("BendingInitialConfiguration","Bending initial configuration M_y",Rfilter, modelMaxMy.analysis, modelMaxMy.halfSegmentNelems, max_segments(5), 1, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("BendingInitialConfiguration","Bending initial configuration M_z",Rfilter, modelMaxMz.analysis, modelMaxMz.halfSegmentNelems, max_segments(6), 1, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("TorsionInitialConfiguration","Torsion initial configuration M_s",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, max_segments(4), 1, cutTreshold, penal, false);
% 

Rfilter = 1.5*(R-r);
penal=3;
cutTreshold = 0.01;
volFr=0.4;

% [xopt_av, xopt_av_lambda, xoptBuckling_av, xoptBuckling_av_lambda, xopt_max, xopt_max_lambda, xoptBuckling_max, xoptBuckling_max_lambda, newLoadFactor ] = configurationTopologyMulti(E,nu,R,r,res, res_thickness, segmentLength,ShapeFn,alpha,modeSamples,frameElems, real_max_segments, loadFactor,Rfilter, cutTreshold, penal, volFr, false, 01);
% 
% save("ComposedTopologyMultiMaxAvOneRing_Rev.mat");
% load("ComposedTopologyMultiMaxAvOneRing_Rev.mat");

% 
% disp(['Average lambda linear = ' num2str(xopt_av_lambda) ' Average  lambda buckling = ' num2str(xoptBuckling_av_lambda)]);
% disp(['Envelope lambda linear   = ' num2str(xopt_max_lambda) ' Envelope lambda buckling   = ' num2str(xoptBuckling_max_lambda)]);
% disp(['New load factor         = ' num2str(newLoadFactor)]);
% 
% developResultsForTopology("OneRingAverageLinear","Average topology linear",modelMaxMs, R, r, alpha, res , segmentLength, segmentNo,xopt_av,sampleMaxMz,dispFactor,Rfilter,cutTreshold,penal);
% developResultsForTopology("OneRingEnvelopeLinear","Envelope topology linear",modelMaxMs, R, r, alpha, res , segmentLength, segmentNo,xopt_max,sampleMaxMz,dispFactor,Rfilter,cutTreshold,penal);
% developResultsForTopology("OneRingAverageBuckling","Average topology with buckling",modelMaxMs, R, r, alpha, res , segmentLength, segmentNo,xoptBuckling_av,sampleMaxMz,dispFactor,Rfilter,cutTreshold,penal);
% developResultsForTopology("OneRingEnvelopeBuckling","Envelope topology with buckling",modelMaxMs, R, r, alpha, res , segmentLength, segmentNo,xoptBuckling_max,sampleMaxMz,dispFactor,Rfilter,cutTreshold,penal);
% 
% plotArmTopOptConfigProjections("OneRingAverage","Aaverage topology",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xoptBuckling_av, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("OneRingEnvelope","Envelope topology",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xoptBuckling_max, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("OneRingAverageLinear","Aaverage topology linear",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xopt_av, cutTreshold, penal, false);
% plotArmTopOptConfigProjections("OneRingEnvelopeLinear","Envelope topology linear",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xopt_max, cutTreshold, penal, false);

[xopt_av, xopt_av_lambda, xoptBuckling_av, xoptBuckling_av_lambda, xopt_max, xopt_max_lambda, xoptBuckling_max, xoptBuckling_max_lambda, newLoadFactor, topOptAvLinear, topOptEnvLinear, topOptAvBuckling, topOptEnvBuckling ] = configurationTopologyMulti(E,nu,R,r,res, res_thickness, segmentLength,ShapeFn,alpha,modeSamples,frameElems, real_max_segments, loadFactor,Pz, Rfilter, cutTreshold, penal, volFr, false, 11);

save("ComposedTopologyMultiMaxAvTwoRingsTh3_Rev.mat","xopt_av", "xopt_av_lambda", "xoptBuckling_av", "xoptBuckling_av_lambda", "xopt_max", "xopt_max_lambda", "xoptBuckling_max", "xoptBuckling_max_lambda", "newLoadFactor", "topOptAvLinear", "topOptEnvLinear", "topOptAvBuckling", "topOptEnvBuckling" );
%load("ComposedTopologyMultiMaxAvTwoRingsTh3_Rev.mat","xopt_av", "xopt_av_lambda", "xoptBuckling_av", "xoptBuckling_av_lambda", "xopt_max", "xopt_max_lambda", "xoptBuckling_max", "xoptBuckling_max_lambda", "newLoadFactor", "topOptAvLinear", "topOptEnvLinear", "topOptAvBuckling", "topOptEnvBuckling");

axisDir   = [0 0 1];
axisPoint = [R 0 0];  % center in XY

%topOptAvLinear.FEAnalysis.mesh

% plotUnwrappedCylMesh(topOptAvLinear.FEAnalysis.mesh.nodes, topOptAvLinear.FEAnalysis.mesh.elems, axisPoint, axisDir, ShapeFn, ...
%     'FieldElem', xopt_av, ...
%     'FaceColor', [0.8 0.8 0.8], ...
%     'FaceAlpha', 1.0);

disp(['Average lambda linear = ' num2str(xopt_av_lambda) ' Average  lambda buckling = ' num2str(xoptBuckling_av_lambda)]);
disp(['Envelope lambda linear   = ' num2str(xopt_max_lambda) ' Envelope lambda buckling   = ' num2str(xoptBuckling_max_lambda)]);
disp(['New load factor         = ' num2str(newLoadFactor)]);

modelMaxN.analysis.felems{1}.edge_color='k';
modelMaxN.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxN.analysis.felems{1}.face_alpha=1;
modelMaxTy.analysis.felems{1}.edge_color='k';
modelMaxTy.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxTy.analysis.felems{1}.face_alpha=1;
modelMaxTz.analysis.felems{1}.edge_color='k';
modelMaxTz.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxTz.analysis.felems{1}.face_alpha=1;
modelMaxMy.analysis.felems{1}.edge_color='k';
modelMaxMy.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxMy.analysis.felems{1}.face_alpha=1;
modelMaxMz.analysis.felems{1}.edge_color='k';
modelMaxMz.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxMz.analysis.felems{1}.face_alpha=1;
modelMaxMs.analysis.felems{1}.edge_color='k';
modelMaxMs.analysis.felems{1}.face_color=[0.8 0.8 0.8];
modelMaxMs.analysis.felems{1}.face_alpha=1;

% numResults=zeros(2,2,4,6);

ConfigNames = {'N';'T_y';'T_z';'M_s';'M_y'; 'M_z'};

maxHM_dd = zeros(6,1);
maxDisp_dd = zeros(6,1);
maxHM_top = zeros(6,1);
maxDisp_top = zeros(6,1);

% plotStep5('segment_model_N', 'segment model, configuration max N',   R, modelMaxN,  topOptEnvLinear, [-1 0 0 0 0 0] );
% plotStep5('segment_model_Ty', 'segment model, configuration max Ty', R, modelMaxTy, topOptEnvLinear, [0 1 0 0 0 0] );
% plotStep5('segment_model_Tz', 'segment model, configuration max Tz', R, modelMaxTz, topOptEnvLinear, [0 0 1 0 0 0] );
% plotStep5('segment_model_Ms', 'segment model, configuration max Ms', R, modelMaxMs, topOptEnvLinear, [0 0 0 1 0 0] );
% plotStep5('segment_model_My', 'segment model, configuration max My', R, modelMaxMy, topOptEnvLinear, [0 0 0 0 1 0] );
% plotStep5('segment_model_Mz', 'segment model, configuration max Mz', R, modelMaxMz, topOptEnvLinear, [0 0 0 0 0 1] );

% plotSegmentResultsStep6('segment_model_AvLin', 'Linear average model', topOptAvLinear, 0.4);
% plotSegmentResultsStep6('segment_model_EnvLin', 'Linear envelope model', topOptEnvLinear, 0.4);
% plotSegmentResultsStep6('segment_model_AvBuckling', 'Linear average model', topOptAvBuckling, 0.4);
% plotSegmentResultsStep6('segment_model_EnvBuckling', 'Linear average model', topOptEnvBuckling, 0.4);

% plotSegmentResultsStep7('segment_topology_AvLin', 'Average linear', topOptAvLinear, newLoadFactor, volFr);
% plotSegmentResultsStep7('segment_topology_EnvLin', 'Envelope linear', topOptEnvLinear, newLoadFactor, volFr );
% plotSegmentResultsStep7('segment_topology_AvBuck', 'Average buckling', topOptAvBuckling, newLoadFactor, volFr);
% plotSegmentResultsStep7('segment_topology_EnvBuck', 'Envelope buckling', topOptEnvBuckling, newLoadFactor, volFr);

% plotSegmentResultsStep8('arm_topology_AvLin', 'Average linear', modelMaxN, topOptAvLinear, volFr);

[maxHM_dd(1), maxDisp_dd(1), maxHM_top(1), maxDisp_top(1)] = prepareResults("linear_average","maxN", modelMaxN, R, r, alpha, res , res_thickness, segmentLength, max_seg_N, xopt_av, topOptAvLinear, sampleMaxN, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(2), maxDisp_dd(2), maxHM_top(2), maxDisp_top(2)] = prepareResults("linear_average","maxT_y", modelMaxTy, R, r, alpha, res , res_thickness, segmentLength, max_seg_Ty, xopt_av, topOptAvLinear, sampleMaxTy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(3), maxDisp_dd(3), maxHM_top(3), maxDisp_top(3)] = prepareResults("linear_average","maxT_z", modelMaxTz, R, r, alpha, res , res_thickness, segmentLength, max_seg_Tz, xopt_av, topOptAvLinear, sampleMaxTz, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(4), maxDisp_dd(4), maxHM_top(4), maxDisp_top(4)] = prepareResults("linear_average","maxM_s", modelMaxMs, R, r, alpha, res , res_thickness, segmentLength, max_seg_Ms, xopt_av, topOptAvLinear, sampleMaxMs, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(5), maxDisp_dd(5), maxHM_top(5), maxDisp_top(5)] = prepareResults("linear_average","maxM_y", modelMaxMy, R, r, alpha, res , res_thickness, segmentLength, max_seg_My, xopt_av, topOptAvLinear, sampleMaxMy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(6), maxDisp_dd(6), maxHM_top(6), maxDisp_top(6)] = prepareResults("linear_average","maxM_z", modelMaxMz, R, r, alpha, res , res_thickness, segmentLength, max_seg_Mz, xopt_av, topOptAvLinear, sampleMaxMz, dispFactor,Rfilter,cutTreshold,penal,Pz);

T1 = table(ConfigNames,maxHM_dd,maxDisp_dd,maxHM_top,maxDisp_top);
T1.Properties.VariableNames = {'Config name' 'max HM design domain'  'u_max design domain'  'max HM topology'  'u_max topology'  };
T1.Properties.Description = 'Linear analysis, average stress intensity';
[maxHM_dd(1), maxDisp_dd(1), maxHM_top(1), maxDisp_top(1)] = prepareResults("linear_envelope","maxN", modelMaxN, R, r, alpha, res , res_thickness, segmentLength, max_seg_N, xopt_max, topOptEnvLinear, sampleMaxN, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(2), maxDisp_dd(2), maxHM_top(2), maxDisp_top(2)] = prepareResults("linear_envelope","maxT_y", modelMaxTy, R, r, alpha, res , res_thickness, segmentLength, max_seg_Ty, xopt_max, topOptEnvLinear, sampleMaxTy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(3), maxDisp_dd(3), maxHM_top(3), maxDisp_top(3)] = prepareResults("linear_envelope","maxT_z", modelMaxTz, R, r, alpha, res , res_thickness, segmentLength, max_seg_Tz, xopt_max, topOptEnvLinear, sampleMaxTz, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(4), maxDisp_dd(4), maxHM_top(4), maxDisp_top(4)] = prepareResults("linear_envelope","maxM_s", modelMaxMs, R, r, alpha, res , res_thickness, segmentLength, max_seg_Ms, xopt_max, topOptEnvLinear, sampleMaxMs, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(5), maxDisp_dd(5), maxHM_top(5), maxDisp_top(5)] = prepareResults("linear_envelope","maxM_y", modelMaxMy, R, r, alpha, res , res_thickness, segmentLength, max_seg_My, xopt_max, topOptEnvLinear, sampleMaxMy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(6), maxDisp_dd(6), maxHM_top(6), maxDisp_top(6)] = prepareResults("linear_envelope","maxM_z", modelMaxMz, R, r, alpha, res , res_thickness, segmentLength, max_seg_Mz, xopt_max, topOptEnvLinear, sampleMaxMz, dispFactor,Rfilter,cutTreshold,penal,Pz);

T2 = table(ConfigNames,maxHM_dd,maxDisp_dd,maxHM_top,maxDisp_top);
T2.Properties.VariableNames = {'Config name' 'max HM design domain'  'u_max design domain'  'max HM topology'  'u_max topology'  };
T2.Properties.Description = 'Linear analysis, envelope stress intensity';

[maxHM_dd(1), maxDisp_dd(1), maxHM_top(1), maxDisp_top(1)] = prepareResults("buckling_average","maxN", modelMaxN, R, r, alpha, res, res_thickness , segmentLength, max_seg_N, xoptBuckling_av, topOptAvBuckling, sampleMaxN, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(2), maxDisp_dd(2), maxHM_top(2), maxDisp_top(2)] = prepareResults("buckling_average","maxT_y", modelMaxTy, R, r, alpha, res , res_thickness, segmentLength, max_seg_Ty, xoptBuckling_av, topOptAvBuckling, sampleMaxTy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(3), maxDisp_dd(3), maxHM_top(3), maxDisp_top(3)] = prepareResults("buckling_average","maxT_z", modelMaxTz, R, r, alpha, res, res_thickness , segmentLength, max_seg_Tz, xoptBuckling_av, topOptAvBuckling, sampleMaxTz, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(4), maxDisp_dd(4), maxHM_top(4), maxDisp_top(4)] = prepareResults("buckling_average","maxM_s", modelMaxMs, R, r, alpha, res, res_thickness , segmentLength, max_seg_Ms, xoptBuckling_av, topOptAvBuckling, sampleMaxMs, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(5), maxDisp_dd(5), maxHM_top(5), maxDisp_top(5)] = prepareResults("buckling_average","maxM_y", modelMaxMy, R, r, alpha, res, res_thickness , segmentLength, max_seg_My, xoptBuckling_av, topOptAvBuckling, sampleMaxMy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(6), maxDisp_dd(6), maxHM_top(6), maxDisp_top(6)] = prepareResults("buckling_average","maxM_z", modelMaxMz, R, r, alpha, res, res_thickness , segmentLength, max_seg_Mz, xoptBuckling_av, topOptAvBuckling, sampleMaxMz, dispFactor,Rfilter,cutTreshold,penal,Pz);

T3 = table(ConfigNames,maxHM_dd,maxDisp_dd,maxHM_top,maxDisp_top);
T3.Properties.VariableNames = {'Config name' 'max HM design domain'  'u_max design domain'  'max HM topology'  'u_max topology'  };
T3.Properties.Description = 'Linear analysis with buckling, average stress intensity';

[maxHM_dd(1), maxDisp_dd(1), maxHM_top(1), maxDisp_top(1)] = prepareResults("buckling_envelope","maxN", modelMaxN, R, r, alpha, res , res_thickness, segmentLength, max_seg_N, xoptBuckling_max, topOptEnvBuckling, sampleMaxN, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(2), maxDisp_dd(2), maxHM_top(2), maxDisp_top(2)] = prepareResults("buckling_envelope","maxT_y", modelMaxTy, R, r, alpha, res , res_thickness, segmentLength, max_seg_Ty, xoptBuckling_max, topOptEnvBuckling, sampleMaxTy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(3), maxDisp_dd(3), maxHM_top(3), maxDisp_top(3)] = prepareResults("buckling_envelope","maxT_z", modelMaxTz, R, r, alpha, res , res_thickness, segmentLength, max_seg_Tz, xoptBuckling_max, topOptEnvBuckling, sampleMaxTz, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(4), maxDisp_dd(4), maxHM_top(4), maxDisp_top(4)] = prepareResults("buckling_envelope","maxM_s", modelMaxMs, R, r, alpha, res , res_thickness, segmentLength, max_seg_Ms, xoptBuckling_max, topOptEnvBuckling, sampleMaxMs, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(5), maxDisp_dd(5), maxHM_top(5), maxDisp_top(5)] = prepareResults("buckling_envelope","maxM_y", modelMaxMy, R, r, alpha, res , res_thickness, segmentLength, max_seg_My, xoptBuckling_max, topOptEnvBuckling, sampleMaxMy, dispFactor,Rfilter,cutTreshold,penal,Pz);
[maxHM_dd(6), maxDisp_dd(6), maxHM_top(6), maxDisp_top(6)] = prepareResults("buckling_envelope","maxM_z", modelMaxMz, R, r, alpha, res , res_thickness, segmentLength, max_seg_Mz, xoptBuckling_max, topOptEnvBuckling, sampleMaxMz, dispFactor,Rfilter,cutTreshold,penal,Pz);

T4 = table(ConfigNames,maxHM_dd,maxDisp_dd,maxHM_top,maxDisp_top);
T4.Properties.VariableNames = {'Config name' 'max HM design domain'  'u_max design domain'  'max HM topology'  'u_max topology'  };
T4.Properties.Description = 'Linear analysis with buckling, envelope stress intensity';

% developResultsForTopology("TwoRingsEnvelope","Envelope topology",modelMaxMs, R, r, alpha, res , segmentLength, segmentNo,xoptBuckling_max,sampleMaxMs,dispFactor,Rfilter,cutTreshold,penal);
% developResultsForTopology("TwoRingsAverageLinear","Average topology linear",modelMaxMs, R, r, alpha, res , segmentLength, segmentNo,xopt_av,sampleMaxMs,dispFactor,Rfilter,cutTreshold,penal);
% developResultjsForTopology("TwoRingsEnvelopeLinear","Envelope topology linear",modelMaxMs, R, r, alpha, res , segmentLength, segmentNo,xopt_max,sampleMaxMs,dispFactor,Rfilter,cutTreshold,penal);
% 
plotArmTopOptConfigProjections("TwoRingsAverage", "Average topology",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xoptBuckling_av, cutTreshold, penal, false);
plotArmTopOptConfigProjections("TwoRingsEnvelope","Envelope topology",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xoptBuckling_max, cutTreshold, penal, false);
plotArmTopOptConfigProjections("TwoRingsAverageLinear", "Average topology linear",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xopt_av, cutTreshold, penal, false);
plotArmTopOptConfigProjections("TwoRingsEnvelopeLinear","Envelope topology linear",Rfilter, modelMaxMs.analysis, modelMaxMs.halfSegmentNelems, 2, xopt_max, cutTreshold, penal, false);


T1.Properties.Description
T1

T2.Properties.Description
T2

T3.Properties.Description
T3

T4.Properties.Description
T4


disp(['Execution finished']);
                                                                                                                                                                                                                                                                      