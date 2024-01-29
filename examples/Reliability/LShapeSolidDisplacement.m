clear;
close all;
l=3;
c=0.4;
res = 14;
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


xp=[l 0.2*l 0.4*l];
P = [2 0 -10];

model = LShapeSolidModel(ShapeFunctionL8,l,res,E,nu,xp,[l 0.2*l 0.4*l]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),0.2) RandomVariable("Normal",P(2),0.1) RandomVariable("Normal",P(3),1)};
transform=IndependentTransformation(randomVariables);
g=loadSolidDisplacementPerformanceFunction(model,fatigueData);

%g.tabNCycles(0.1,500,10)

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.15;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;

tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000000, 2);
tuner.tuneMC();
tuner.plotMCs(["Px" "Py" "Pz"],'Nc');

%tuner.tuneFORM();

 sora2 = SORA('LShapeSolidDispBeta_2', model,topOpt, randomVariables, g, transform, 2);
 sora3 = SORA('LShapeSolidDispBeta_3', model,topOpt, randomVariables, g, transform, 3);
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

