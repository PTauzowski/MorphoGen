clear;
close all;
l=3;
c=0.4;
res = 20;
E=210000;
nu=0.3;
betat=3;

xp=[l 0.4*l];
P = [0 -10];

model = LShapeModelLinear(ShapeFunctionL4,l,res,E,nu,xp,P);
model.setResultNode([0 l]);
model.plotModel();

randomVariables={RandomVariable("Normal",P(1),5) RandomVariable("Normal",P(2),5)};
transform=IndependentTransformation(randomVariables);
g=loadPenalizedStressPerformanceFunction(model);
g.computeValue(P)
hmv = HMV(randomVariables,g,transform,betat);
form = FORM(randomVariables,g,transform);

[ Pf_form, mpp_form, beta_form ]= form.solve();
[ Pf_hmv, mpp_hmv, beta_hmv ] = hmv.solve();
fprintf("\nDesign Domain \n BetaDet_form=%5.7f,  BetaDet_hmv=%5.7f\n",beta_form,beta_hmv);
if beta_form > betat
    fprintf("\nDesign domain reliability indes lower than target beta\n");
    return;
end

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.4;

topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
sora = SORA(model, topOpt, randomVariables, g, transform, betat );


Pf_mc=0;
% N=5000;
% mc= MonteCarlo(randomVariables,g,N);
% [ Pf_mc, p ] = mc.solve();
% mc.scatterPlots(["Px" "Py"],"HMpen");

[ Pf_form, mpp_form, beta_form ]= form.solve();
[ Pf_hmv, mpp_hmv, beta_hmv ] = hmv.solve();
fprintf("\nDeterministic design \n Pf_Det_mc=%5.7f, BetaDet_form=%5.7f,  BetaDet_hmv=%5.7f\n",Pf_mc,beta_form,beta_hmv);
if beta_form > betat
    fprintf("\nDesign domain reliability lower than target beta\n");
end

g.computeValue(P)
sora.limitReliability();

 % [betaDet, xRel, volDet, volRel] = sora.solve();
 % figure;
 % plot(sora.mpps(:,1),sora.mpps(:,2),'Marker','o');
 % [ Pf, mpp, betar ] = form.solve()
