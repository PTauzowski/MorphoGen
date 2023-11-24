clear;
close all;
res = 20;
l = 3;
sf = ShapeFunctionL4;
endAngle = 45;
E=1;
nu=0.3;
P=[-1 0];

model = InclinedSupportPlanestress(sf,l,endAngle,res,E,nu,P);

model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.0);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode([l 0]);
model.computeDisplacement(1,0.3,[0 -1])

model.setResultNode([0 0]);
model.computeHMstress(1,0.3,[0 -1])



