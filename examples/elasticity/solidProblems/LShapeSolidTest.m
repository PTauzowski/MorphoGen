clear;
close all;
res = 10;
E=210000;
nu=0.3;
l = 3;
xp=[l 0.2*l 0.4*l];
P=[0 0 -100];

model = LShapeSolidModel(ShapeFunctionL8,l,res,E,nu,xp,P);

model.plotModel();
model.solveWeighted();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode([l 0 0]);
model.computeDisplacement(E,nu,P)
model.setResultNode([l 0.4*l 0]);
model.computeDisplacement(E,nu,P)
model.setResultNode([l 0 0]);
model.computeHMstress(E,nu,P)
model.setResultNode([l 0.4*l 0]);
model.computeHMstress(E,nu,P)
model.setResultNode([0.4*l 0 0.4*l]);
model.computeHMstress(E,nu,P)
model.setResultNode([0.4*l 0.4*l 0.4*l]);
model.computeHMstress(E,nu,P)



