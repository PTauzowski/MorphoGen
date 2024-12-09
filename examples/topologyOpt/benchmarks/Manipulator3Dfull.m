clear;
close all;

E=210E09;
nu=0.3;
R=0.25;
r=0.22;
segmentLength=0.3;
res=15;
alpha=30;

nArms=7;
nSamples=5000;
samples=random("Uniform",0,360,nSamples,nArms);
samples(:,1)=0;
maxNM=zeros(nSamples,1);    
endPoints=zeros(nSamples,3);
ShapeFn = ShapeFunctionL8;

% tic;
% for k=1:nSamples
%     betas=samples(k,:);
% 
%     model = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, betas, ShapeFn);
%     model.analysis.printProblemInfo();
% 
%     x=ones(model.analysis.getTotalElemsNumber(),1);
%     model.analysis.solveWeighted(x);
%     model.analysis.computeElementResults();
%     maxHM(k)= max(model.fe.results.nodal.all(:,13));
%     endPoints(k,:)=model.xEnd;
% end
% disp(['Average single analysis time: ' num2str(toc/nSamples)]);

load("ManipulatorOpti5000.mat");

r=0.23;
res=10;

[vMin, imin]=min(maxHM);
[vMax, imax]=max(maxHM);
[vSort, isort]=sort(maxHM);

modelMin = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(imin,:), ShapeFn);
% modelMin.fe.plot(modelMin.mesh.nodes);
% modelMin.analysis.plotCurrentLoad();
%line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
% title('Model for minimal Huber-Mises');

% figure
% x=ones(modelMin.analysis.getTotalElemsNumber(),1);
% modelMin.analysis.solveWeighted(x);
% modelMin.analysis.computeElementResults();
% modelMin.analysis.plotMaps(["sHM"],0.1);
% modelMin.fe.plotWired(modelMin.mesh.nodes,modelMin.analysis.qnodal,0.1);

figure;
modelMax = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(isort(nSamples-2),:), ShapeFn);
%modelMax.fe.plot(modelMax.mesh.nodes);
% modelMin.analysis.plotCurrentLoad();
%line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
% title('Model for maximal Huber-Mises');

% figure
% x=ones(modelMax.analysis.getTotalElemsNumber(),1);
% modelMax.analysis.solveWeighted(x);
% modelMax.analysis.computeElementResults();
% modelMax.analysis.plotMaps(["sHM"],0.1);
% modelMax.fe.plotWired(modelMax.mesh.nodes,modelMax.analysis.qnodal,0.1);



% nSort=12;
% 
% for k=1:nSort
%     sortId = nSamples-k+1;
%     %sortId = nSamples/nSort*k;
%     modelSort = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, samples(isort(sortId),:), ShapeFn);
%     %modelSort.fe.plot(modelMin.mesh.nodes);
%     %line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
% 
%     x=ones(modelSort.analysis.getTotalElemsNumber(),1);
%     modelSort.analysis.solveWeighted(x);
%     modelSort.analysis.computeElementResults();
%     modelSort.analysis.plotMaps(["sHM"],0.1);
%     modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
%     title(['Model for minimal Huber-Mises for HMmax=' num2str(vSort(sortId))]);
% 
% end


% frameElems=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
% mesh=Mesh();
% frameElem=Frame3D(frameElems,E,0.02,0.8*E,0.0004,0.0004,0.003);
% 
% figure;
% frameElem.plot(modelMin.frameNodes);
% mesh.nodes=modelMin.frameNodes;
% 
% analysis = LinearElasticityWeighted( frameElem, mesh, false );
% analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 0 0 0 0]);
% 
% analysis.loadClosestNode(mesh.nodes(end,:), ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 -1 0 0 0] );
% analysis.plotCurrentLoad();
% analysis.plotSupport();
% 
% x=ones(model.analysis.getTotalElemsNumber(),1);
% analysis.solveWeighted(x);
% [Fel, Feg] = frameElem.computeResults(mesh.nodes,analysis.qnodal);
% 
% Fel
% 
% % Row and column descriptions
% rowNames = {
%     'Axial Force Start', 'Shear Force Y Start', 'Shear Force Z Start', ...
%     'Torsion Moment X Start', 'Moment Y Start', 'Torsion Moment Z Start', ...
%     'Axial Force End', 'Shear Force Y End', 'Shear Force Z End', ...
%     'Torsion Moment X End', 'Moment Y End', 'Moment Z End'
% };
% 
% columnNames = {'Bar 1', 'Bar 2', 'Bar 3', 'Bar 4', 'Bar 5', 'Bar 6', 'Bar 7'};
% 
% % Create table
% internalForcesTable = array2table(Fel, 'VariableNames', columnNames, 'RowNames', rowNames);
% % Add title as a property
% internalForcesTable.Properties.Description = 'Internal Forces for 3D Frame Finite Elements';
% 
% % Display the table
% disp(internalForcesTable);
 
% Filtering radius
Rfilter = 1.5*(R-r);
penal=3;
cutTreshold = 0.02;

% tic
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, modelMax.analysis, cutTreshold, penal, 0.4, false );
% topOpt.setConstElems(modelMax.const_elems);
% [objF, xopt]  = topOpt.solve();
% toc
% 
% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, modelMax.analysis, penal, 0.3, false);
% [objF, xopt]  = topOpt.solve();
% toc
% 
% 


