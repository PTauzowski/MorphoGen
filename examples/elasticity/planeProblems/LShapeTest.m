clear;
close all;

% four noded Lagrange shape functions
sf=ShapeFunctionL4;

% Length of the longest edge
l=1;

% mesh resolution along shortest edge
nl=60;

%Young modulus
E=1;

%Poisson's ratio
nu=0.3;

% Loacation of a force
xp=[l 0.4*l];

% Value of the acting force 
P=[0 -1];

% LShape model object instantiation
model = LShapeModelLinear(sf,l,nl,E,nu,xp,P);

% solving linear elastic problem
model.solveWeighted();

%plotting of contour maps
model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

%computing displacements for selected node
model.setResultNode([0 0]);
model.computeDisplacement(1,0.3,[0 -1])

%computing Huber-Mises stress for selected node
model.setResultNode([0.4*l 0.4*l]);
model.computeHMstress(1,0.3,[0 -1])
