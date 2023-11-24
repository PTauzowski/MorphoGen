clear;
close all;

sf=ShapeFunctionL4;
l=1;
nl=60;
E=1;
nu=0.3;
xp=[l 0.4*l];
P=[0 -1];

model = LShapeModelLinear(sf,l,nl,E,nu,xp,P);

model.solve();

model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode([0 0]);
model.computeDisplacement(1,0.3,[0 -1])

model.setResultNode([0.4*l 0.4*l]);
model.computeHMstress(1,0.3,[0 -1])
