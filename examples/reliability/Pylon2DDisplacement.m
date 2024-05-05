clear;
close all;
h=15;
h1=3;
b=1.5;
b1=10;
c=0.4;
resb = 30;
E=210000;
nu=0.3;

xp=[0 h-h1];
P = [0 -10];

model = Pylon2DModel(ShapeFunctionL4,h,h1,b,b1,resb,E,nu,xp,[-b/2 h-h1]);
model.plotModel();
% model.solveWeighted();
% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
% model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

randomVariables={RandomVariable("Normal",P(1),0.2) RandomVariable("Normal",P(2),0.1) };
transform=IndependentTransformation(randomVariables);
g = performanceFunctionPlus( displacementPerformanceFunction(model,2) );

%g.tabNCycles(0.1,500,10)

Rfilter = 1.2*b/resb;
penal = 3;
cutTreshold = 0.001;
volFr=0.4;
% 
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, false );
topOpt.is_silent=true;
topOpt.solve();

% tic
% topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, model.analysis, penal, volFr, false);
% [objF, xopt]  = topOpt.solve();
%  toc

% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py" "Pz"],'Nc');

%tuner.tuneFORM();
g.threshold=0.0;

 sora2 = SORA('PylonPlaneDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('PylonPlaneDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);
 %sora.checkTuning();

 
% topOpt.solve();
%sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
% sora2.tabReliability();

%sora2.checkTuning();

% sora_results2 = sora2.solveX();
% sora_results3 = sora3.solveX();

% form_res = form.solve()

