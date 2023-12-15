clear;
close all;
randomVariables={RandomVariable("Normal",0.4,0.1) RandomVariable("Normal",0.4,0.1)};
transform=IndependentTransformation(randomVariables);
g=SinCosApr();
N=100000;
mc = MonteCarlo(randomVariables,g,N);
form = FORM(randomVariables,g,transform);
mc_results = mc.solve();
fprintf( "\nMonte Carlo results: Pf=%1.6f, beta=%1.6f\n",mc_results.Pf, norminv( mc_results.Pf ) );
form_results = form.solve();
fprintf( "FORM results: Pf=%1.6f, beta=%1.6f\n",form_results.Pf, form_results.beta );

mc.scatterPlots(["x1" "x2"],"objFn");

