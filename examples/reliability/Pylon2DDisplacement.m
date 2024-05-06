clear;
close all;
h=15;
h1=3;
b=1.5;
b1=10;
c=0.4;
resb = 15;
E=210000;
nu=0.3;

xp=[b1 h-h1];
xres=[0 h-h1];
P = [0 -10];

model = Pylon2DModel(ShapeFunctionL4,h,h1,b,b1,resb,E,nu,xp,xres);
model.plotModel();
% model.solveWeighted();
% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
% model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

randomVariables={RandomVariable("Normal",P(1),0.5) RandomVariable("Normal",P(2),0.5)  };
transform=IndependentTransformation(randomVariables);
g = performanceFunctionPlus( displacementPerformanceFunction(model,2) );

%g.tabNCycles(0.1,500,10)

Rfilter = 1.2*b/resb;
penal = 3;
cutTreshold = 0.001;
volFr=0.3;
% 
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, false );
topOpt.is_silent=true;

% model.setupLoad([P(1) P(2) P(1) P(2)]);
% topOpt.solve();
model.plotModel();

% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, model.analysis, penal, volFr, false);
% [objF, xopt]  = topOpt.solve();
%  toc


 tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 100000000, 2);
 % tuner.checkTopology(Pdest);
 % [topMin,topMax] = tuner.getBoundsTopologies(0.3,0.6);
 % tuner.fullReliabilityTuning(Pdest,topMin,topMax);

g.threshold = 0.0025;

 % 
 % tuner.tuneMC();
 % tuner.plotMCs(["Px" "Py"],'u');


 %mpps = tuner.checkModality(0.5)
%tuner.tuneFORM();


 sora2 = SORA('PylonPlaneDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('PylonPlaneDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);
 sora2 = SORA('PylonPlaneDispBeta20_4', model,topOpt, randomVariables, g, transform, 4);
 sora3 = SORA('PylonPlaneDispBeta20_5', model,topOpt, randomVariables, g, transform, 5);
 %sora.checkTuning();

 
% topOpt.solve();
%sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
% sora2.tabReliability();

%sora2.checkTuning();

 sora_results2 = sora2.solveX();
% sora_results3 = sora3.solveX();

% form_res = form.solve()

