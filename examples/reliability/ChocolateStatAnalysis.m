clear;
close all;
height=12;
E=210000;
nu=0.3;
alphaT=1 ;
dT=-20;

 % x(1) - gan, alGan proportion,  x(2) - relNotchDepth ,    x(3) - notchWidth
 % 0.2358    0.4383    5.0757 - test
lb=[0.2 0.1 1];
ub=[0.8 0.9 8];
x0=(lb+ub)/2;

m=(ub+lb)/2;
s=sqrt(((ub-lb).^2/12));

randomVariables={RandomVariable("Uniform",lb(1),ub(1)) RandomVariable("Uniform",lb(2),ub(2)) RandomVariable("Uniform",lb(3),ub(3))};

transform=IndependentTransformation(randomVariables);
g = chocolatePerformanceFunction(height,210000,0.3,alphaT,dT);
fn_g = @(x)( g.computeValue(x) );
%g.fullFactorialBoundsPlot(lb,ub);
N=100;
mc= MonteCarlo(randomVariables,g,N);
%x = mc.generateRandomSapmles(N);
tic
%res_mc = mc.solve();
xopt = fmincon(fn_g,x0,[],[],[],[],lb,ub)
%xopt = fmincon(fn_g,x0)
toc

%save("chocolateStat5000_dT.mat");

%[v, i]=max(mc.r)
%g.evaluateValue(mc.x(i,:))
g.evaluateValue(xopt);
%g.createModel([0.2358    0.4383    5.0757]);
%g.evaluateValue2();



