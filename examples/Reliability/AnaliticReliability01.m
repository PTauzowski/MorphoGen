clear;
close all;
randomVariables={RandomVariable("Normal",0.4,0.1) RandomVariable("Normal",0.4,0.1)};
g=SinCosApr();
N=1000000;
mc = MonteCarlo(randomVariables,g,N);
form = FORM(randomVariables,g);
[Pfmc p] = mc.solve();
Pfmc
Pform = form.solve()

figure, hold on;
scatter3(mc.x(p>0,1),mc.x(p>0,2),p(p>0),'MarkerEdgeColor',[0 .8 .8],'Marker','.');
scatter3(mc.x(p<=0,1),mc.x(p<=0,2),p(p<=0),'filled','MarkerEdgeColor',[0.5 0 .5],'Marker','o');
