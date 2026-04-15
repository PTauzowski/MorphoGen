clear; clc;
close all;

% ======================================================================
%  PILLAR MODEL — PARAMETER STRUCT
%  All lengths in the same unit (e.g. µm or mm — be consistent).
%  The model exploits quarter symmetry: the XZ- and YZ-planes are
%  treated as symmetry planes, so only one quarter of the full geometry
%  is meshed.
% ======================================================================

% ---- Shape function --------------------------------------------------
p.sf = ShapeFunctionH27();
%  H27: 27-node tri-quadratic hexahedral element (2nd-order accuracy).
%  Replace with ShapeFunctionH8() for a 1st-order (8-node) mesh.

% ---- Pillar geometry -------------------------------------------------
p.top_R = 50;
%  Radius of the pillar at its top surface.

p.pillar_inclination_deg = 5;
%  Taper half-angle measured from the vertical axis (degrees).
%  0 = straight cylinder.  5 gives a gentle inward taper toward the top
%  (Rbase = top_R + pillar_height * tan(inclination_deg)).

p.pillar_layers = [ 10   400   30   2.6   30 ];
%  Layer thicknesses of the pillar, ordered BOTTOM → TOP.
%  Layer  1 ( 10):  lower cap / mounting section
%  Layer  2 (400):  main pillar body
%  Layer  3 ( 30):  upper transition
%  Layer  4 (2.6):  thin functional layer (e.g. piezoelectric)
%  Layer  5 ( 30):  top cap

p.pillar_res = [ 1   8   2   1   2 ];
%  Number of elements in the z-direction for each pillar layer
%  (same order as pillar_layers).

% ---- Ground / substrate geometry ------------------------------------
p.ground_layers = [ 200   50   300 ];
%  Layer thicknesses of the ground substrate, ordered BOTTOM → TOP.
%  Layer 1 (200):  deep substrate
%  Layer 2 ( 50):  transition zone
%  Layer 3 (300):  top ground layer (contains the depression ring)

p.ground_res = [ 10   2   10 ];
%  Number of elements in the z-direction for each ground layer
%  (same order as ground_layers).

% ---- Interface (pillar–ground boundary) -----------------------------
p.int_th = 0.1;
%  Thickness of the transitional interface slab that straddles z = 0
%  (before the z_offset shift is applied).  The slab spans
%  z ∈ [−int_th/2, +int_th/2] and takes int_th/2 from the top ground
%  layer and int_th/2 from the bottom pillar layer.

% ---- Material chemistry (per layer) ---------------------------------
p.pillar_chem = [ 0.00   0.08   0.08   0.18   0.08 ];
%  Pillar chemistry value assigned to each pillar layer (same order as
%  pillar_layers).  In the current FEAP export this is used as the
%  pillar-side composition amplitude, so the 10 nm GaN cap is 0.00.

p.ground_chem = [ 0   1   1 ];
%  Ground chemistry value assigned to each ground layer (same order as
%  ground_layers).  In the current FEAP export this populates the second
%  chemistry channel used on the ground side.

% ---- Vertical offset ------------------------------------------------
p.z_offset = -p.pillar_layers(1);
%  Rigid shift applied to ALL z-coordinates after mesh generation.
%  Setting it to −pillar_layers(1) places z = 0 at the top of the
%  first pillar layer, i.e. at the base of the main pillar body.

% ---- Depression (moat around the pillar base) -----------------------
%  A smooth bowl-shaped depression is applied to the top surface of the
%  inner annular ground ring surrounding the pillar base.  Leave a field
%  as NaN (or omit it) to let the constructor choose a default that is
%  derived from the geometry.

p.depression_width = NaN;
%  Radial width of the depression ring  (Rout − Rin).
%  NaN → 20 % of Rbase (the pillar base radius).

p.depression_depth = NaN;
%  Maximum depth of the depression at its lowest point.
%  NaN → 4 % of the effective top ground-layer thickness, then
%        safety-capped so no element inversion occurs.

p.depression_r_min = NaN;
%  Radial coordinate r at which the depression reaches its minimum.
%  NaN → inner edge + 25 % of depression_width.
%  Must satisfy  Rin < depression_r_min < Rout.

% ---- Tile (outer square boundary) -----------------------------------
%  Both pillar diameters (100 nm and 2000 nm) share the same array pitch
%  (same lithography grid).  Tile size is therefore fixed by the wide-pillar
%  geometry (top_R = 1000) and reused for the small pillar.
%  Formula mirrors the NaN auto rule: 1.2 × Rbank, Rbank = 1.5 × Rout,
%  Rout ≈ top_R_wide + pillar_height × tan(inclination).
pillar_height_total  = sum( p.pillar_layers );
Rout_wide            = 1000 + pillar_height_total * tand( p.pillar_inclination_deg );
p.tile_size          = 1.2 * 1.5 * Rout_wide;
%  Half-width of the outer square tile (distance from the symmetry axis
%  to the outer boundary along x or y).
%  Must be strictly larger than Rbank; otherwise the constructor errors.

% ---- Mesh resolution ------------------------------------------------
p.res_cyl  = 4;
%  Circumferential node divisions per 90° arc used for every circular
%  cross-section (pillar cylinder, ground annuli, bank rings).
%  Increase for a finer angular discretisation.

p.res_ring = 5;
%  Radial node divisions of each annular (ring / pipe) section.
%  Controls the element count from inner to outer radius of each ring.

p.res_tile = 20;
%  Radial node divisions of the square tile transition sections.
%  Increase independently of res_ring for a finer outer region.

% ======================================================================
%  BUILD MODELS
% ======================================================================

model = PillarModel( p );

% ---- Wide-tile variant -----------------------------------------------
%  Same pillar geometry, embedded in a much larger substrate tile.
%  All pillar / ground / chemistry / interface parameters are inherited
%  from p; only the tile size and tile resolution differ.
p_wide           = p;
p_wide.top_R     = 1000;   % wider pillar top radius for the second model
p_wide.tile_size = NaN;    % auto: 1.2 × Rbank (larger tile, consistent proportions)
p_wide.res_tile  = 15;     % finer tile resolution for the larger domain

p_wide.depression_width = 120;
%  Radial width of the depression ring  (Rout − Rin).
%  NaN → 20 % of Rbase (the pillar base radius).

p_wide.res_ring = 10;
%  Radial node divisions of each annular (ring / pipe) section.
%  Controls the element count from inner to outer radius of each ring.

model_width = PillarModel( p_wide );

% ======================================================================
%  DIAGNOSTICS
% ======================================================================

[badE, minDetJ] = model.mesh.findNegativeJacobian(p.sf, 1e-12);
fprintf("Bad elements before: %d (min detJ = %.3e)\n", numel(badE), min(minDetJ));

% ======================================================================
%  VISUALISATION
% ======================================================================

fe  = SolidElasticElem( p.sf, model.mesh.elems );
fe2 = SolidElasticElem( p.sf, model_width.mesh.elems );

model.plotLayerColors(fe);

figure;
model_width.plotLayerColors(fe2);

% ======================================================================
%  EXPORT
% ======================================================================

model.FEAP_Export("Pillar18.i");
model_width.FEAP_Export("Pillar18w.i");
