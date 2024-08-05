clear;
close all;
sf=ShapeFunctionL4;
l=1;
nl=60;
E=1;
nu=0.3;
xp=[l l/4];
P=[0 -1];

model = CantileverModelLinear(sf,l,nl,E,nu,xp);
model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

