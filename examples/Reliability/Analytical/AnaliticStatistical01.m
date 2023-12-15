clear;
close all;
randomVariables={RandomVariable("Normal",0.7,0.25) RandomVariable("Normal",0.7,0.25)};
g=SinCosG();
N=100000;
stat = StatisticalAnalysis(randomVariables,g);
r = stat.solve(N);
scatter(r.x(:,1),r.x(:,2));


