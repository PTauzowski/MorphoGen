clear;
close all;
l=3;
c=0.4;
res = 40;
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
fatigueData.Nexp=42150;


xp=[l 0.4*l];
P = [0 10];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp,P);
model.setResultNode([l*0.4 l*0.4]);
model.plotModel();

randomVariables={RandomVariable("Normal",0,2) RandomVariable("Normal",-8,0.2)};
transform=IndependentTransformation(randomVariables);
g=loadFatiguePerformanceFunction(model,fatigueData);

% N=5000;
% mc= MonteCarlo(randomVariables,g,N);
% [ Pf_mc, p ] = mc.solve();
% Pf_mc
% mc.scatterPlots(["Px" "Py"],"Ux");
% % 
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

