clear;
close all;
sf=ShapeFunctionL16;
l=0.6;
h=5;
nl=100;
E=80000;
nu=0.3;
P=[0 -1];
xp=[0 h];

model = ColumnModel(sf,l,h,nl,E,nu,P,xp);
model.analysis.solve(5);

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
figure;
model.analysis.setForm(1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.2);




