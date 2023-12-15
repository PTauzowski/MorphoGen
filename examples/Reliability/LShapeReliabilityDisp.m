clear;
close all;
l=3;
c=0.4;
res = 20;
E=210000;
nu=0.3;

xp=[l 0.4*l];
P = [0 -100];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp,P);
model.setResultNode([l*0.4 l*0.4]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),10) RandomVariable("Normal",P(2),10)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunctionDisp(model);
betat=3;

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.4;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
topOpt.is_silent=true;

tuner = ReliabilityTaskTuner(model, topOpt, randomVariables, transform, g, 5000, betat);

tuner.tuneMC();
%tuner.plotMCs(["Px" "Py"],'Ux');

% 
% sora = SORA(model,topOpt, randomVariables, g, transform, 3);
% sora.limitReliability();
% sora.tabReliability();
% sora_results = sora.solve();
% form_res = form.solve()
% res_mc = mc.solve()
% [xDet, betaDet, xRel, volDet, volRel]  = sora.solve();
% volDet
% volRel


