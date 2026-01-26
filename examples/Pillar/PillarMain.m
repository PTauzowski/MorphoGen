clear; clc;
close all;

top_R = 100;
pillar_layers = [ 10 400 30 2 6 30 ];
pillar_res    = [ 1  8  2  1 1 2 ];
ground_layers = [ 200 50 300 ];
ground_res    = [ 10   2 10 ];
pillar_chem = [0.08 0.08 0.08 0.08 0.18  0.08 ];
ground_chem = [0 1 1];
int_th = 0.1;

sf = ShapeFunctionH27();

model = PillarModel( top_R, pillar_layers, ground_layers, pillar_res, ground_res, pillar_chem, ground_chem, int_th, sf );
model.checkMeshIntegrity(0.001);

% diagnose
[badE, minDetJ] = model.mesh.findNegativeJacobian(sf, 1e-12);
fprintf("Bad elements before: %d (min detJ = %.3e)\n", numel(badE), min(minDetJ));

% fix by renumbering
res = model.mesh.fixNegativeJacobianByRenumbering(sf, 1e-12);
disp(res)

% re-check
[badE2, minDetJ2] = model.mesh.findNegativeJacobian(sf, 1e-12);
fprintf("Bad elements after: %d (min detJ = %.3e)\n", numel(badE2), min(minDetJ2));
model.checkMeshIntegrity(0.001);

figure;
fe = SolidElasticElem( sf, model.mesh.elems );
fe.plotWithSettings(model.mesh.nodes);

figure;
detJ = fe.computeDets( model.mesh.nodes ,[] );
bad_elems=find(detJ(1,1,:,14)<0);
feb = SolidElasticElem( sf, model.mesh.elems(bad_elems,:) );
feb.plot(model.mesh.nodes, [0.8,0,0]);

filename = "Pillar11.i";
model.FEAP_Export(filename);

[nodes, elements] = readFEAPFile(filename);
fe1 = SolidElasticElem( sf, elements );
figure;
fe1.plot(nodes);