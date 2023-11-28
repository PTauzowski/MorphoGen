clear;
close all;

ganTh=10;
alGanTh=5;
notchWidth=4;
relNotchDepth=0.875;
relRoutndNotchDepth=0.3;

E=210000;
nu=0.3;
pressure=100;

model = ChocolateModel( ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth, E, nu);

%model.plotModel();
model.solveWeighted();
%model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
%model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.computeStressObjective()

