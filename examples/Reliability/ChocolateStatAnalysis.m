clear;
close all;
height=15;
E=210000;
nu=0.3;

randomVariables={RandomVariable("Uniform",0.2,0.8) RandomVariable("Uniform",0.2,0.8), RandomVariable("Uniform",2,8)};
transform=IndependentTransformation(randomVariables);
g = chocolatePerformanceFunction(height,210000,0.3);

 N=1000;
 mc= MonteCarlo(randomVariables,g,N);
 tic
 res_mc = mc.solve();
 toc
 res_mc.Pf

 save("chocolateStat1000_2.mat");




