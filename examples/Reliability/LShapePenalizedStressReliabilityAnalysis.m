clear;
close all;
l=3;
c=0.4;
res = 20;
E=210000;
nu=0.3;

xp=[l 0.4*l];
P = [0 10];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp,P);
model.setResultNode([0 l]);
model.plotModel();

randomVariables={RandomVariable("Normal",-10,1) RandomVariable("Normal",-10,1)};
transform=IndependentTransformation(randomVariables);
g=loadPenalizedStressPerformanceFunction(model);

% N=1000;
% mc= MonteCarlo(randomVariables,g,N);
% [ Pf_mc, p ] = mc.solve();
% Pf_mc
% mc.scatterPlots(["Px" "Py"],"Ux");
% 
 hmv = HMV(randomVariables,g,transform,3);
% form = FORM(randomVariables,g,transform);
% Pf_form = form.solve()
% [ Pf, mpp, betar ] = hmv.solve()

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.4;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
%[objF, xopt]  = topOpt.solve();

sora = SORA(model,topOpt, hmv);
[xDet, betaDet, xRel, volDet, volRel] = sora.solve();
volDet
volRel

