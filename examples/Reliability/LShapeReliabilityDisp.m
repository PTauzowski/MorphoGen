clear;
close all;
l=3;
c=0.4;
res = 80;
E=210000;
nu=0.3;

xp=[l 0.4*l];
P = [0 -10];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp);
model.setResultNode([l l*0.4]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),3) RandomVariable("Normal",P(2),1)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionDisp(model);
%gl=loadLinearPerformanceFunction(model,g);

% g.computeValue([1 0])
% g.computeValue([0 1])
% 
% gl.computeValue([1 0])
% gl.computeValue([0 1])


Rfilter = l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.25;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
%topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, model.analysis, penal, volFr, true);
[objF, xopt]  = topOpt.solve();
toc

topOpt.is_silent=true;

% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000000, 2);
% tic;
% tuner.tuneMC();
% toc
% tuner.plotMCs(["Px" "Py"],'HMstress');
%tuner.tuneFORM();

 sora2 = SORA('LShapeDispBeta20_2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('LShapeDispBeta20_3', model,topOpt, randomVariables, g, transform, 3);
 sora4 = SORA('LShapeDispBeta20_4', model,topOpt, randomVariables, g, transform, 4);
 sora5 = SORA('LShapeDispBeta20_5', model,topOpt, randomVariables, g, transform, 5);

 % sora2.solveX();
 % sora3.solveX();
 %  sora4.solveX();
 % sora5.solveX();
 
 % topOpt.solve();
%  sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
% sora2.tabReliability();
% 
% sora2.checkTuning();
% sora3.checkTuning();
   


