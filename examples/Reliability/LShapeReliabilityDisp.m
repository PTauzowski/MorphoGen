clear;
close all;
l=3;
c=0.4;
res = 40;
E=210000;
nu=0.3;

xp=[l 0.4*l];
P = [0 -100];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp);
model.setResultNode([l l*0.4]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),10) RandomVariable("Normal",P(2),10)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionDisp(model);
%gl=loadLinearPerformanceFunction(model,g);

% g.computeValue([1 0])
% g.computeValue([0 1])
% 
% gl.computeValue([1 0])
% gl.computeValue([0 1])


Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.25;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;

% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000, 2);
% tic;
% tuner.tuneMC();
% toc
% tuner.plotMCs(["Px" "Py"],'Ux');
%tuner.tuneFORM();

 sora2 = SORA('LShapeDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('LShapeDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);

 
 % topOpt.solve();
%  sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
%  sora2.tabReliability();

%sora2.checkTuning();
   
sora_results2 = sora2.solveX();
sora_results3 = sora3.solveX();

