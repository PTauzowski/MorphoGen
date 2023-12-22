clear;
close all;
l=3;
c=0.4;
res = 40;
E=210000;
nu=0.3;

xp=[l 0.4*l];
P = [0 -100];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp,P);
model.setResultNode([l l*0.4]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),10) RandomVariable("Normal",P(2),10)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionDisp(model);

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.001;
volFr=0.3;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;

tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 5000, 2);
%tuner.tuneMC();
%tuner.plotMCs(["Px" "Py"],'Ux');
%tuner.tuneFORM();

 sora2 = SORA('LShapeDispBeta2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('LShapeDispBeta3', model,topOpt, randomVariables, g, transform, 3);
 %sora.checkTuning();

 
 % topOpt.solve();
%  sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
 %sora.tabReliability();

 %sora2.checkTuning();

  sora_results2 = sora2.solveX();
  sora_results3 = sora3.solveX();
% form_res = form.solve()

