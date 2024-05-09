clear;
close all;
res = 100;
l = 1;

aspect=2;

Rfilter = 3*l/res;

%Removal intensity threshold
cutTreshold = 0.005;

%penalty factor
penal = 3;

sfL4 = ShapeFunctionL4;
mesh = Mesh();
mesh.addRectMesh2D(0, 0, aspect*l, l, aspect*res, res, sfL4.pattern);
fe=PlaneStressElem( sfL4, mesh.elems );

material = PlaneStressMaterial('mat1');
material.setElasticIzo(1, 0.3);
fe.setMaterial( material );
fe.props.h=1;

analysis = LinearElasticityWeighted( fe, mesh, true );
fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.001 ) );
loadedEdgeSelector = Selector( @(x)( (abs( x(:,1) - aspect*l) < 1.0E-4) & ( abs(x(:,2)-l)<l/15 )  ) );
analysis.loadClosestNode([aspect*l, l/2 ], ["ux" "uy"], [0 -1] );
analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );

analysis.printProblemInfo();
fe.plotSolid(mesh.nodes);
analysis.plotCurrentLoad();
analysis.plotSupport();


% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, 0.05, penal, 0.4, true );
% tic
% [objF, xopt]  = topOpt.solve();
% title(['cutThreshold=' num2str(topOpt.maxais) ', ' num2str(topOpt.iteration) ' iterations, execution time ' num2str(toc)]);
% 
% figure;
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, 0.04, penal, 0.4, true );
% tic
% [objF, xopt]  = topOpt.solve();
% title(['cutThreshold=' num2str(topOpt.maxais) ', '  num2str(topOpt.iteration) ' iterations, execution time ' num2str(toc)]);
% 
% figure;
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, 0.03, penal, 0.4, true );
% tic
% [objF, xopt]  = topOpt.solve();
% title(['cutThreshold=' num2str(topOpt.maxais) ', ' num2str(topOpt.iteration) ' iterations, execution time ' num2str(toc)]);
% 
% figure;
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, 0.02, penal, 0.4, true );
% tic
% [objF, xopt]  = topOpt.solve();
% title(['cutThreshold=' num2str(topOpt.maxais) ', ' num2str(topOpt.iteration) ' iterations, execution time ' num2str(toc)]);
% 
% figure;
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, 0.01, penal, 0.4, true );
% tic
% [objF, xopt]  = topOpt.solve();
% title(['cutThreshold=' num2str(topOpt.maxais) ', ' num2str(topOpt.iteration) ' iterations, execution time ' num2str(toc)]);
% 
% figure;
% topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, 0.005, penal, 0.4, true );
% tic
% [objF, xopt]  = topOpt.solve();
% title(['cutThreshold=' num2str(topOpt.maxais) ', ' num2str(topOpt.iteration) ' iterations, execution time ' num2str(toc)]);

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, 0.05, penal, 0.4, true );
topOpt.ct1=0.07;
topOpt.ct2=0.001;
tic
[objF, xopt]  = topOpt.solve();
title(['cutThreshold variable =' num2str(topOpt.ct1) '-' num2str(topOpt.ct2) ', ' num2str(topOpt.iteration) ' iterations, execution time ' num2str(toc)]);

% figure;
% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, analysis, penal, 0.4, true);
% [objF, xopt]  = topOpt.solve();
% toc



