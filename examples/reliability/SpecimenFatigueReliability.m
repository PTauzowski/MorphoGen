clear;
close all;
a=10;
div=15;
c=0.4;
E=200000;
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
%fatigueData.Fref=25;
fatigueData.Nexp=42150;

%P = 60;
P = 150;

model = SpecimenModelLinear(ShapeFunctionL4,a,div,E,nu,P);
% model.setResultNode([0.4*l 0.2*l 0.4*l]);
model.plotModel();

randomVariables={RandomVariable("Normal",E,0.15*E) RandomVariable("Normal",fatigueData.sy,0.15*fatigueData.sy)};
transform=IndependentTransformation(randomVariables);
g=SpecimenFatiguePerformanceFunction(model,fatigueData);


%g.tabNCycles(50,350,100)

Rfilter = 3*30/1000/div/3;
penal = 3;
cutTreshold = 0.005;
volFr=0.75;

%topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance(Rfilter, model.analysis, penal, volFr, false);
topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, false );
topOpt.is_silent=true;
topOpt.const_elems=g.model.const_elems;
model.setX(topOpt.x);

 
% tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 1000, 2);
% tuner.tuneMC();
% tuner.plotMCs(["Px" "Py" "Pz"],'Nc');

%tuner.tuneFORM();

 sora2 = SORA('SpecimenFatigueBeta4', model,topOpt, randomVariables, g, transform, 4);
 sora3 = SORA('SpecimenFatigueBeta5', model,topOpt, randomVariables, g, transform, 5);
 %sora.checkTuning();

 
 %topOpt.solve();
% sora2.limitReliability();
 %sora.tabMultiMpp();

% topOpt.solve();
% sora2.tabReliability();

%sora2.checkTuning();

sora_results2 = sora2.solveX();
sora_results3 = sora3.solveX();

% form_res = form.solve()

