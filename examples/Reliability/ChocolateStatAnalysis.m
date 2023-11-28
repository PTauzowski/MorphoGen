clear;
close all;
height=15;
E=210000;
nu=0.3;

randomVariables={RandomVariable("Uniform",0.2,0.8) RandomVariable("Uniform",0.2,0.8), RandomVariable("Uniform",2,8)};
transform=IndependentTransformation(randomVariables);
g = chocolatePerformanceFunction(height,210000,0.3);

 N=10;
 mc= MonteCarlo(randomVariables,g,N);
 [ Pf_mc, p ] = mc.solve();
 Pf_mc
 figure, hold on;
 scatter3(mc.x(p>0,1),mc.x(p>0,2),p(p>0),'MarkerEdgeColor',[0 .8 .8],'Marker','.');
 scatter3(mc.x(p<=0,1),mc.x(p<=0,2),p(p<=0),'filled','MarkerEdgeColor',[0.5 0 .5],'Marker','o');




