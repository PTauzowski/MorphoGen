clear;
close all;
r1=5;
r2=10;
div=20;
sf = ShapeFunctionT3;
E=210000;
nu=0.3;
P=[0 -1];


model = LameProblemModelTriangular(sf,r1,r2,div,E,nu,P);

model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.0);
%model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode([r1 0]);
model.computeDisplacement(1,0.3,[0 -1])
model.setResultNode([r2 0]);
model.computeDisplacement(1,0.3,[0 -1])

model.setResultNode([r1 0]);
model.computeHMstress(1,0.3,[0 -1])
model.setResultNode([r2 0]);
model.computeHMstress(1,0.3,[0 -1])
