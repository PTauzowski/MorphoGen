clear;
close all;
w=0.5;
l=8*w;
h=4*w;
c=0.4;
res = 10;
E=210000;
nu=0.3;

fatigueData.E = E;
fatigueData.sy = 190;
fatigueData.Cy = 6000;
fatigueData.su = 450;
fatigueData.epD = 0.12;
fatigueData.m=2;
fatigueData.sfi = 140;
fatigueData.S = 2.8;
fatigueData.s = 2;
fatigueData.Dc = 0.2;
fatigueData.Fref=0.01;
fatigueData.Nexp=42150;


xp=[l w/2 0];
P = [-2 0 -5];


model = CantileverSolidModelLinear(ShapeFunctionL8,w,res,E,nu,xp,xp);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),0.2) RandomVariable("Normal",P(2),0.1) RandomVariable("Normal",P(3),0.5)};
transform=IndependentTransformation(randomVariables);
g=loadCantileverSolidPerformanceFunctionDisp(model);

%g.tabNCycles(0.1,500,10)

Rfilter = 2*w/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.1;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;

% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py" "Pz"],'Ux');

%tuner.tuneFORM();

 sora2 = SORAold('CantileverSolidDispBeta_2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORAold('CantileverSolidDispBeta_3', model,topOpt, randomVariables, g, transform, 3);
 %sora.checkTuning();

 
% topOpt.solve();
%sora2.limitReliability();
  %sora.tabMultiMpp();

% topOpt.solve();
% sora2.tabReliability();

 % sora2.checkTuning();
 % sora3.checkTuning();

 sora_results2 = sora2.solveX();
 sora_results3 = sora3.solveX();

% form_res = form.solve()

