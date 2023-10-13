randomVariables={RandomVariable("Normal",0.7,0.25) RandomVariable("Normal",0.7,0.25)};
g=SinCosApr();
N=100000;
mc = MonteCarlo(g,randomVariables,N);
form = FORM(g,randomVariables);
[Pfmc p] = mc.solve();
Pfmc
Pform = form.solve()

figure, hold on;
scatter3(mc.x(p>0,1),mc.x(p>0,2),p(p>0),'MarkerEdgeColor',[0 .8 .8],'Marker','.');
scatter3(mc.x(p<=0,1),mc.x(p<=0,2),p(p<=0),'filled','MarkerEdgeColor',[0.5 0 .5],'Marker','o');

