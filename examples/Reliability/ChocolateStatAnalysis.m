clear;
close all;
height=15;
E=210000;
nu=0.3;
alphaT=11E-1;
dT=20;
lb=[0.3 0.6 2];
ub=[0.9 0.8 8];

m=(ub+lb)/2;
s=sqrt(((ub-lb).^2/12));

randomVariables={RandomVariable("Uniform",lb(1),ub(1)) RandomVariable("Uniform",lb(2),ub(2)) RandomVariable("Uniform",lb(3),ub(3))};

transform=IndependentTransformation(randomVariables);
g = chocolatePerformanceFunction(height,210000,0.3,alphaT,dT);
%g.fullFactorialBoundsPlot(lb,ub);
N=5000;
mc= MonteCarlo(randomVariables,g,N);
%x = mc.generateRandomSapmles(N);
tic
res_mc = mc.solve();
toc

save("chocolateStat200_dT.mat");

[v, i]=max(mc.r)
%g.evaluateValue(mc.x(i,:));
g.evaluateValue([0.67,0.7,2]);



