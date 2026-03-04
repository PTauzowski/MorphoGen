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

model = PillarModel( top_R, pillar_layers, ground_layers, pillar_res, ground_res, pillar_chem, ground_chem, int_th, sf , z_offset);
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

figure;
fe = SolidElasticElem( sf, model.mesh.elems );

% --- layer-coloured plot ---
% Element z-centroids (vectorised, works for H27 with 27 nodes).
% Linear indexing v(matrix) preserves matrix shape; 2-D indexing A(matrix,col)
% does not — MATLAB linearises the row-subscript, giving a column vector.
node_z = model.mesh.nodes(:, 3);                        % 47439 × 1
cz     = mean(node_z(model.mesh.elems), 2);             % nElems × 1

% Canonical layer-boundary z-planes shifted to exported coordinate system
zP = computeCanonicalZ(int_th, ground_layers, pillar_layers) + z_offset;

% Colour palette — blues for ground layers, warm hues for pillar layers.
% Add rows here if you add more layers.
layer_colors = [
    0.15, 0.35, 0.75;   % ground 1  – dark blue
    0.20, 0.55, 0.90;   % ground 2  – medium blue
    0.30, 0.70, 0.95;   % ground 3  – light blue
    0.95, 0.55, 0.10;   % pillar 1  – orange
    0.85, 0.15, 0.15;   % pillar 2  – red
    0.15, 0.70, 0.25;   % pillar 3  – green
    0.90, 0.80, 0.10;   % pillar 4  – yellow
    0.55, 0.10, 0.80;   % pillar 5  – purple
];
iface_color = [0.78 0.78 0.78];   % grey for interface slabs

hold on;  daspect([1 1 1]);  axis off;  view(3);

nominal_k = 0;
for k = 1 : numel(zP) - 1
    z0   = zP(k);
    z1   = zP(k+1);
    dz   = z1 - z0;
    mask = cz >= z0 - 1e-9 & cz <= z1 + 1e-9;
    if ~any(mask), continue; end

    if abs(dz - int_th) < int_th * 0.1   % interface slab → grey
        col = iface_color;
    else                                   % nominal layer  → colour
        nominal_k = nominal_k + 1;
        col = layer_colors(min(nominal_k, size(layer_colors,1)), :);
    end

    fe.plotWithSettings(model.mesh.nodes, ...
        "elem nums", find(mask), "elem color", col, "edge color", "k");
end

% figure;
% detJ = fe.computeDets( model.mesh.nodes ,[] );
% bad_elems=find(detJ(1,1,:,14)<0);
% feb = SolidElasticElem( sf, model.mesh.elems(bad_elems,:) );
% feb.plot(model.mesh.nodes, [0.8,0,0]);

filename = "Pillar16.i";
model.FEAP_Export(filename);

% Reading and displaying mesh from FEAP file
% [nodes, elements] = readFEAPFile(filename);
% fe1 = SolidElasticElem( sf, elements );
% figure;
% fe1.plot(nodes);

% =========================================================================
function zPlanes = computeCanonicalZ(int_th, ground_layers, pillar_layers)
% Canonical element-face z-planes for a PillarModel with z_offset = 0.
% Mirrors the h_eff arithmetic of addLayeredQuarterCylinder / addLayeredPipe3D.
%
%   Ground stack : z = -sum(ground_layers) … -int_th/2
%   Global seam  : z = -int_th/2 … +int_th/2
%   Pillar stack : z = +int_th/2 … sum(pillar_layers)

    Lg = numel(ground_layers);
    Lp = numel(pillar_layers);

    % Ground effective inputs (last layer trimmed for global seam)
    ge     = ground_layers(:);
    ge(Lg) = ge(Lg) - int_th/2;

    % Internal h_eff for ground (same rule as addLayeredPipe3D)
    hg = ge;
    for k = 1:Lg
        if k > 1,  hg(k) = hg(k) - int_th/2; end
        if k < Lg, hg(k) = hg(k) - int_th/2; end
    end

    % Pillar effective inputs (first layer trimmed for global seam)
    pe    = pillar_layers(:);
    pe(1) = pe(1) - int_th/2;

    % Internal h_eff for pillar
    hp = pe;
    for k = 1:Lp
        if k > 1,  hp(k) = hp(k) - int_th/2; end
        if k < Lp, hp(k) = hp(k) - int_th/2; end
    end

    % Accumulate z-planes bottom → top
    z      = -sum(ground_layers);
    zPlanes = z;

    for k = 1:Lg
        z = z + hg(k);
        zPlanes(end+1) = z; %#ok<AGROW>
        if k < Lg
            z = z + int_th;       % internal ground interface
            zPlanes(end+1) = z;   %#ok<AGROW>
        end
    end
    % z == -int_th/2 here; add global seam top
    z = z + int_th;
    zPlanes(end+1) = z;           %#ok<AGROW>

    for k = 1:Lp
        z = z + hp(k);
        zPlanes(end+1) = z;       %#ok<AGROW>
        if k < Lp
            z = z + int_th;       % internal pillar interface
            zPlanes(end+1) = z;   %#ok<AGROW>
        end
    end
    % z == sum(pillar_layers) here

    zPlanes = sort(unique(zPlanes));
end