clear;
close all;
l=3;
c=0.4;
res = 16;
E=210000;
nu=0.3;

xp=[l 0.2*l 0.4*l];
P = [0 0 -100];

model = LShapeSolidModel(ShapeFunctionL8,l,res,E,nu,xp,P);
model.setResultNode([l 0.2*l 0.4*l]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),10) RandomVariable("Normal",P(2),1) RandomVariable("Normal",P(3),10)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionSolidDisp(model);

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.1;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;
 % topOpt.solve();

% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 5000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py" "Pz"],'Uz');
%tuner.tuneFORM();

 sora2 = SORA('LShapeSolidDispBeta2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('LShapeSolidDispBeta3', model,topOpt, randomVariables, g, transform, 3);
 %sora.checkTuning();
 
  %topOpt.solve();
%  sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
 %sora.tabReliability();

 sora2.checkTuning();

  % sora_results2 = sora2.solveX();
  % sora_results3 = sora3.solveX();
% form_res = form.solve()

