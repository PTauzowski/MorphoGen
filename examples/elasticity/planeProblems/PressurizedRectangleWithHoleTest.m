clear;
close all;
res = 20;
l = 3;
x0 = [5 5];
a=10;
hf = 0.5;
sf = ShapeFunctionL4;
E=210000;
nu=0.3;
pressure=100;

model = PressurizedSquareWithHoleModel(sf,x0,a,hf,res,E,nu,pressure);

model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.0);
%model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode(x0-[a/2*hf 0]);
model.computeDisplacement(E,nu,pressure)
model.setResultNode(x0+[a/2*hf 0]);
model.computeDisplacement(E,0.3,pressure)

model.setResultNode(x0-[a/2*hf 0]);
model.computeHMstress(E,0.3,pressure)
model.setResultNode(x0+[a/2*hf 0]);
model.computeHMstress(E,0.3,pressure)
