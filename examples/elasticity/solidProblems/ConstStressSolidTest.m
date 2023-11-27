clear;
close all;

res = 20;
l = 3;
E=210000;
nu=0.3;
pressure=100;

modelL8  = ConstStressSolidModel(ShapeFunctionL8,l,res,E,nu,pressure);
modelL27 = ConstStressSolidModel(ShapeFunctionL27,l,res,E,nu,pressure);
%modelL64 = ConstStressSolidModel(ShapeFunctionL64,l,res,E,nu,pressure);
%modelT6  = ConstStressSolidModel(ShapeFunctionT6,l,res,E,nu,pressure);

model=modelL8;

model.plotModel();
model.solveWeighted();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode([l 0 0]);
model.computeDisplacement(E,nu,pressure)
model.setResultNode([l l/4 0]);
model.computeHMstress(E,nu,pressure)
model.setResultNode([l/2 l/4 l/4]);
model.computeHMstress(E,nu,pressure)

