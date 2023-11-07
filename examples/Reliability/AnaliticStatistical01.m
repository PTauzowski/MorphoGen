clear;
close all;
randomVariables={RandomVariable("Normal",0.7,0.25) RandomVariable("Normal",0.7,0.25)};
g=SinCosG();
N=100000;
stat = StatisticalAnalysis(randomVariables,g);
sc = stat.solve(N);
scatter(sc(:,1),sc(:,2));

% theta = linspace(0,1,500);
% x = exp(theta).*sin(100*theta);
% y = exp(theta).*cos(100*theta);
% s = scatter(x,y);


