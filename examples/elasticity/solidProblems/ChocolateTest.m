clear;
close all;

ganTh=4;
alGanTh=8;
notchWidth=4;
relNotchDepth=0.875;
relRoutndNotchDepth=0.3;

E=210000;
nu=0.3;
alphaT=11.0E-6;
dT = 20;
pressure=100;

model = ChocolateModel( ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth, E, nu, alphaT, dT);

model.plotModel();
model.plotZCoordsPoints()

% model.solveWeighted();
% model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
% model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);
% 
% model.computeStressObjective()

