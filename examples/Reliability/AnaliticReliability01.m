randomVariables={RandomVariable("Normal",0.7,0.25) RandomVariable("Normal",0.7,0.25)};
g=SinCosG();
N=100000;
mc = MonteCarlo(g,randomVariables,N);
form = FORM(g,randomVariables);
Pfmc = mc.solve()
Pform = mc.solve()

