clear;
close all;

E=210E09;
nu=0.3;
R=0.25;
r=0.20;
segmentLength=0.3;
res=10;
alpha=30;

nArms=6;
nSamples=5000;
samples=random("Uniform",0,360,nSamples,nArms);
maxNM=zeros(nSamples,1);    
endPoints=zeros(nSamples,3);
ShapeFn = ShapeFunctionL8;

tic;
for k=1:nSamples
    betas=samples(k,:);
   
    model = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, betas, ShapeFn);
    model.analysis.printProblemInfo();
    
    x=ones(model.analysis.getTotalElemsNumber(),1);
    model.analysis.solveWeighted(x);
    model.analysis.computeElementResults();
    maxHM(k)= max(model.fe.results.nodal.all(:,13));
    endPoints(k,:)=model.xEnd;
end
disp(['Average single analysis time: ' num2str(toc/nSamples)]);

[vMin, imin]=min(maxHM);
[vMax, imax]=max(maxHM);

modelMin = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(imin,:), ShapeFn);
modelMin.fe.plot(modelMin.mesh.nodes);
line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
title('Model for minimal Huber-Mises');

figure;
modelMax = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(imax,:), ShapeFn);
modelMax.fe.plot(modelMax.mesh.nodes);
line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
title('Model for maximal Huber-Mises');

%save("ManipulatorOpti5000.mat");

%problem.plotNodes();
%model.analysis.plotCurrentLoad();
%model.analysis.plotSupport();

% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
% model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

 
% Filtering radius
Rfilter = 1.5*(R-r);
penal=3;
cutTreshold = 0.05;

% tic
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, true );
% [objF, xopt]  = topOpt.solve();
% toc
% 
% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, true);
% [objF, xopt]  = topOpt.solve();
% toc
% 
% 


