clear;
close all;
height=12;
E=210000;
nu=0.3;
alphaT=11E-6;
dT=-20;
lb=[0.3 0.6 3];
ub=[0.9 0.8 8];
x0=(lb+ub)/2;

m=(ub+lb)/2;
s=sqrt(((ub-lb).^2/12));

randomVariables={RandomVariable("Uniform",lb(1),ub(1)) RandomVariable("Uniform",lb(2),ub(2)) RandomVariable("Uniform",lb(3),ub(3))};

transform=IndependentTransformation(randomVariables);
g = chocolatePerformanceFunction(height,210000,0.3,alphaT,dT);
fn_g = @(x)( g.computeValue(x) );
%g.fullFactorialBoundsPlot(lb,ub);
N=1000;
%mc= MonteCarlo(randomVariables,g,N);
%x = mc.generateRandomSapmles(N);
tic
%res_mc = mc.solve();
%xopt = fmincon(fn_g,x0,[],[],[],[],lb,ub);
toc

%save("chocolateStat5000_dT.mat");

%[v, i]=max(mc.r)
%g.evaluateValue(mc.x(i,:))
g.createModel(x0);
%g.evaluateValue2();



