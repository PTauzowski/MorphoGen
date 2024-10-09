clear;
close all;
sf=ShapeFunctionL4;
l=1;
nl=60;
E=1;
nu=0.3;
xp=[l l/4];
P=[0 -1];

% Define nodes and edges
s = [1 1 2 3];
t = [2 3 4 4];
weights = [5 3 2 4]; % Edge weights

% Create graph with weights
G = graph(s, t, weights);

% Add node names
G.Nodes.Name = {'A', 'B', 'C', 'D'}';

% Plot with labels and weights
p = plot(G, 'EdgeLabel', G.Edges.Weight, 'NodeLabel', G.Nodes.Name);
nodeDegree = degree(G);
disp(nodeDegree);

adjMatrix = adjacency(G);
disp(adjMatrix);

p = plot(G);
layout(p, 'force'); % Force-directed layout

model = CantileverModelLinear(sf,l,nl,E,nu,xp);
model.solveWeighted();

model.plotModel();
model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
model.fe.plotWired(model.mesh.nodes,model.analysis.qnodal,0.1);

model.setResultNode([l 0]);
model.computeDisplacement(1,0.3,[0 -1])

model.setResultNode([0 0]);
model.computeHMstress(1,0.3,[0 -1])


