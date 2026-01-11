clear; clc;
close all;

top_R = 100;
pillar_layers = [ 10 400 30 2 6 30 ];
pillar_res    = [ 1  20  3  2 2 3 ];
ground_layers = [ 200 50 300 ];
ground_res    = [ 20   5  20 ];
pillar_chem = [0.08 0.08 0.08 0.18  0.08 ];
ground_chem = [0 1 1];
int_th = 0.1;

sf = ShapeFunctionH27();

model = PillarModel( top_R, pillar_layers, ground_layers, pillar_res, ground_res, pillar_chem, ground_chem, int_th, sf );

figure;
fe = SolidElasticElem( sf, model.mesh.elems );
fe.plot(model.mesh.nodes);