clear;
close all;

res = 20;
l = 3;
E=210000;
nu=0.3;
pressure=[-100 0];

modelL4  = ConstPlaneStressModel(ShapeFunctionL4,l,res,E,nu,pressure);
modelL9  = ConstPlaneStressModel(ShapeFunctionL9,l,res,E,nu,pressure);
modelL16 = ConstPlaneStressModel(ShapeFunctionL16,l,res,E,nu,pressure);
modelT3 = ConstPlaneStressModelTriangular(ShapeFunctionT3,l,res,E,nu,pressure);

model=modelT3;

model.plotModel();
model.solveWeighted();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode([l 0]);
model.computeDisplacement(E,nu,pressure)
model.setResultNode([l l/4]);
model.computeHMstress(E,nu,pressure)
model.setResultNode([l/2 l/4]);
model.computeHMstress(E,nu,pressure)

