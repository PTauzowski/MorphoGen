clear; clc; close all;

addpath(genpath('../..'));

filename = 'Pillar20.i';   % change to Pillar21w.i for the wide variant

[nodes, elems] = readFEAPFile(filename);

% Filter out Winkler Q9 surface elements — they are padded with trailing
% zeros to fit the 27-column layout, so any row with zeros is not a volume elem.
elems = elems(all(elems > 0, 2), :);

mesh      = Mesh();
mesh.nodes = nodes;
mesh.elems = elems;

sf = ShapeFunctionH27();
fe = SolidElasticElem(sf, mesh.elems);

figure;
fe.plotWithSettings(mesh.nodes, ...
    "elem color", [0.65 0.82 0.95], ...
    "edge color", "k",...
    "nodes",true);
axis off;
view(3);
title(filename, 'Interpreter', 'none');