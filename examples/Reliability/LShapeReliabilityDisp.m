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
g=loadPerformanceFunction(model);

% N=5000;
% mc= MonteCarlo(randomVariables,g,N);
% res_mc = mc.solve()
% mc.scatterPlots(["Px" "Py"],"Ux");
% 
% hmv = HMV(randomVariables,g,transform,3);
% form = FORM(randomVariables,g,transform);
% form_res = form.solve()
% hmv_res = hmv.solve()

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.4;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
%[objF, xopt]  = topOpt.solve();

sora = SORA(model,topOpt, randomVariables, g, transform, 3);
sora.limitReliability();
%sora.tabReliability();
% sora_results = sora.solve();
% form_res = form.solve()
% res_mc = mc.solve()
% [xDet, betaDet, xRel, volDet, volRel]  = sora.solve();
% volDet
% volRel


