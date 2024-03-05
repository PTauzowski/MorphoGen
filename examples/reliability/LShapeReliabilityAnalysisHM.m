clear;
close all;
l=3;
c=0.4;
res = 40;
E=210000;
nu=0.3;

xp=[l 0.4*l];
P = [0 10];

baseName='LShapeRelHM';

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp);
model.setResultNode([0.4*l 0.4*l]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),4) RandomVariable("Normal",P(2),2)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionHM(model);
betat=2;

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.3;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;

% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py"],'HM');

%tuner.tuneFORM();

 sora2 = SORA('LShapeHMBeta2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('LShapeHMBeta3', model,topOpt, randomVariables, g, transform, 3);
 %sora.checkTuning();

 
 % topOpt.solve();
%  sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
% sora2.tabReliability();

 %sora2.checkTuning();

sora_results2 = sora2.solveX();
sora_results3 = sora3.solveX();

% form_res = form.solve()

