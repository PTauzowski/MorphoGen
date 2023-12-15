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
%model.plotModel();

randomVariables={RandomVariable("Normal",P(1),5) RandomVariable("Normal",P(2),5)};
transform=IndependentTransformation(randomVariables);
g=loadPenalizedStressPerformanceFunction(model);
g.computeValue(P)
hmv = HMV(randomVariables,g,transform,betat);
form = FORM(randomVariables,g,transform);

% res_form = form.solve();
% res_hmv = hmv.solve();
% fprintf("\nDesign Domain \n BetaDet_form=%5.7f,  BetaPred_hmv=%5.7f\n",res_form.beta, res_hmv.beta_pred);
% if res_form.beta > betat
%     fprintf("\nDesign domain reliability indes lower than target beta\n");
% end

model.setupLoad(P);

Rfilter = 1.2*l/res;
penal = 3;
cutTreshold = 0.005;
volFr=0.4;

% Pf_mc=0;
% N=5000;
% mc= MonteCarlo(randomVariables,g,N);
% res_mc = mc.solve();
% 
% fprintf("\n minDesign :%1.5f, maxDesign :%1.5f\n",min(mc.r),max(mc.r));
 topOpt = StressIntensityTopologyOptimizationVol( Rfilter, model.analysis, cutTreshold, penal, volFr, true );
 topOpt.is_silent=true;
% model.setupLoad(P);
% topOpt.solve();
sora = SORA(model, topOpt, randomVariables, g, transform, betat );

model.setupLoad(P);
% res_mc = mc.solve();
% mc.scatterPlots(["Px" "Py"],"HMpen");
% fprintf("\n minTop :%1.5f, maxTop :%1.5f\n",min(mc.r),max(mc.r));

% res_form = form.solve();
% res_hmv = hmv.solve();
% fprintf("\nDeterministic design \n Beta_Det_mc=%5.7f, BetaDet_form=%5.7f,  BetaDet_hmv=%5.7f\n",-norminv(res_mc.Pf), res_form.beta, res_hmv.beta_pred);
% if res_form.beta> betat
%     fprintf("\nDesign domain reliability lower than target beta\n");
% end
% 
g.computeValue(P)
sora.limitReliability();

sora_res = sora.solve()
 % figure;
 % plot(sora.mpps(:,1),sora.mpps(:,2),'Marker','o');
 % [ Pf, mpp, betar ] = form.solve()
