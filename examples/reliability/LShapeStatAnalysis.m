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
model.setResultNode([l 0]);
model.plotModel();

randomVariables={RandomVariable("Normal",100,20) RandomVariable("Normal",100,10)};
transform=IndependentTransformation(randomVariables);
g=loadPerformanceFunction(model);

N=1000;
mc= MonteCarlo(randomVariables,g,N);
[ Pf_mc, p ] = mc.solve();
Pf_mc

mc.scatterPlots(["Px" "Py"],"Ux");

hmv = HMV(randomVariables,g,transform,3);
form = FORM(randomVariables,g,transform);
Pf_form = form.solve()
[ Pf, mpp, betar ] = hmv.solve()


