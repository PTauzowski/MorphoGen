clear;
close all;
l=3;
c=0.4;
res = 40;
E=210000;
nu=0.3;
betat=3;

xp=[l 0.4*l];
P = [5 -10];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp,P);
model.setResultNode([0 l]);
%model.plotModel();

randomVariables={RandomVariable("Normal",P(1),2) RandomVariable("Normal",P(2),2)};
transform=IndependentTransformation(randomVariables);
g=loadPenalizedStressPerformanceFunction(model);
Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.25;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;

% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 5000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py"],'Ps');

%tuner.tuneFORM();

 sora2 = SORA('LShapePenalizedStressBeta2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('LShapePenalizedStressBeta3', model,topOpt, randomVariables, g, transform, 3);
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
