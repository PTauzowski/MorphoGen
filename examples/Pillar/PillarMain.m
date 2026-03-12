clear; clc;
close all;

top_R = 100;
pillar_layers = [ 10 400 30 2.6 30 ];
pillar_res    = [ 1  8  2  1 2 ];
ground_layers = [ 200 50 300 ];
ground_res    = [ 10   2 10 ];
pillar_chem = [ 0.08 0.08 0.08 0.18  0.08 ];
ground_chem = [0 1 1];
int_th = 0.1;

z_offset = - pillar_layers(1);

sf = ShapeFunctionH27();

model = PillarModel( top_R, pillar_layers, ground_layers, 1.2, pillar_res, ground_res, 1, pillar_chem, ground_chem, int_th, sf , z_offset);
model_width = PillarModel( top_R, pillar_layers, ground_layers, 3, pillar_res, ground_res, 3, pillar_chem, ground_chem, int_th, sf , z_offset);
%model.checkMeshIntegrity(1e-6);

% diagnose
[badE, minDetJ] = model.mesh.findNegativeJacobian(sf, 1e-12);
fprintf("Bad elements before: %d (min detJ = %.3e)\n", numel(badE), min(minDetJ));

% % fix by renumbering
% res = model.mesh.fixNegativeJacobianByRenumbering(sf, 1e-12);
% disp(res)
% 
% % re-check
% [badE2, minDetJ2] = model.mesh.findNegativeJacobian(sf, 1e-12);
% fprintf("Bad elements after: %d (min detJ = %.3e)\n", numel(badE2), min(minDetJ2));
% model.checkMeshIntegrity(0.001);

fe = SolidElasticElem( sf, model.mesh.elems );

fe2 = SolidElasticElem( sf, model_width.mesh.elems );

% --- layer-coloured plot ---
model.plotLayerColors(fe);

figure;
% --- layer-coloured plot ---
model_width.plotLayerColors(fe2);

% figure;
% detJ = fe.computeDets( model.mesh.nodes ,[] );
% bad_elems=find(detJ(1,1,:,14)<0);
% feb = SolidElasticElem( sf, model.mesh.elems(bad_elems,:) );
% feb.plot(model.mesh.nodes, [0.8,0,0]);

filename = "Pillar15.i";
model.FEAP_Export(filename);

filename = "Pillar15w.i";
model_width.FEAP_Export(filename);

% Reading and displaying mesh from FEAP file
% [nodes, elements] = readFEAPFile(filename);
% fe1 = SolidElasticElem( sf, elements );
% figure;
% fe1.plot(nodes);

