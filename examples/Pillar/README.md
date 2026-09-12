# Pillar

Mesh generator and FE model for a quarter axisymmetric semiconductor nano-pillar.

| | |
|---|---|
| **Branch of record** | `Vibrations` (tagged `paper/vibrations`) |
| **Publication** | *TODO — add the DOI once published* |

This study shares the `Vibrations` branch with the beam-vibration work in
`examples/Vibrations/`; they are two separate studies, which is why they have
two folders.

---

## Pillar Mesh Generator

This directory contains a mesh generator for a quarter axisymmetric pillar-on-ground structure, suitable for finite element analysis of semiconductor nano-pillars.

## Files

- `generatePillarQuarterMesh.m` - Main mesh generator function
- `testPillarMesh.m` - Demo and test script
- `Pillar.m` - Existing class file (separate from mesh generator)
- `README.md` - This file

## Quick Start

```matlab
% Generate mesh with default parameters
mesh = generatePillarQuarterMesh();

% Generate with custom resolution
opts.target_dz = 15;  % Finer vertical resolution
opts.nr_core = 12;     % More radial divisions in core
mesh = generatePillarQuarterMesh(opts);

% Use in FEM analysis
sf = ShapeFunctionH27();
fe = SolidElasticElem(sf, mesh.elems);
% ... continue with analysis setup
```

## Physical Structure

The mesh represents a **200 nm diameter pillar** (at top) on a substrate with the following layer stack (from bottom to top):

### Below z=0:
- 50 nm GaN (non-porous substrate)
- 300 nm porous GaN:Si (ends exactly at z=0)

### Above z=0 (pillar only):
- 10 nm GaN
- 400 nm In₀.₀₈Ga₀.₉₂N
- 30 nm In₀.₀₈Ga₀.₉₂N
- 2.6 nm In₀.₁₈Ga₀.₈₂N (quantum well)
- 30 nm In₀.₀₈Ga₀.₉₂N (top layer)

**Pillar geometry:**
- Cone with 11° inclination from vertical
- 200 nm diameter at top (radius = 100 nm)
- Expands downward: `r_pillar(z) = rTop + (zTop - z) * tan(11°)`

**Bank region:**
- Material exists outside pillar up to height `zBank` (default 90 nm)
- Above `zBank`, only pillar exists
- Smooth trench profile below z=0

## Key Design Principles

### 1. Horizontal Layer Interfaces
**CRITICAL:** All layer interfaces are horizontal planes at fixed z-levels. The fabrication process deposits flat layers, then etching creates the pillar. Etching only modifies the solid/void boundary, never warps layer interfaces.

**Automatic Layer Splitting:** If a physical layer would straddle z=0 or zBank, it is automatically split into segments at these critical boundaries. Split segments maintain the same material ID and proportional mesh resolution.

### 2. Coordinate System
- **Quarter model**: x ≥ 0, y ≥ 0 (exploits symmetry)
- **z=0**: Exactly the top of the porous layer
- **r = sqrt(x² + y²)**: Radial distance

### 3. Mesh Strategy
The domain is decomposed into three radial zones:

- **Core (0 → Rp_zone)**: Cylindrical mesh around pillar
- **Trench (Rp_zone → Rtrench)**: Annular region with etched boundary
- **Outer (Rtrench → Rout)**: Circle-to-square transition to tile boundary

### 4. Interface Sublayers
Between every pair of adjacent physical layers, an interface sublayer is inserted:
- Thickness: `tInterface` (default 0.1 nm)
- Always 1 element thick (`nz = 1`)
- Allows modeling of interface effects

## Parameters

### Pillar Geometry
```matlab
opts.rTop = 100;        % Pillar top radius [nm]
opts.alpha_deg = 11;    % Cone angle from vertical [degrees]
opts.zBank = 90;        % Bank height above z=0 [nm]
```

### Domain Size
```matlab
opts.Rout = 1000;       % Outer domain size [nm]
opts.zSubExtra = 0;     % Extra substrate below stack [nm]
```

### Trench/Bank Profile
```matlab
opts.Rtrench = 300;        % Trench zone outer radius [nm]
opts.trenchWidth = 150;    % Radial width of trench [nm]
opts.trenchDepth = 200;    % Trench depth below z=0 [nm]
opts.Hcurve = 50;          % Vertical smoothing height [nm]
opts.trenchSharpness = 2.0; % Profile sharpness
```

### Interface Sublayers
```matlab
opts.tInterface = 0.1;  % Interface thickness [nm]
```

### Mesh Resolution
```matlab
opts.target_dz = 20;    % Target vertical element size [nm]
opts.nr_core = 8;       % Radial divisions in core
opts.nr_trench = 8;     % Radial divisions in trench
opts.nr_outer = 6;      % Radial divisions in outer zone
opts.ntheta = 8;        % Angular divisions in quarter circle
```

### Element Type
```matlab
opts.sfName = 'H27';    % 'H27' (27-node) or 'H8' (8-node)
```

### Etched Region Mode
```matlab
opts.etch_mode = 'void';  % 'void' (no elements) or 'air_material'
```

## Material IDs

The mesh includes material ID assignments in `mesh.matID`:

| Material ID | Description |
|------------|-------------|
| 1 | Extra substrate (if `zSubExtra > 0`) |
| 2 | GaN substrate (50 nm) |
| 3 | Porous GaN:Si (300 nm) |
| 4 | GaN (10 nm, pillar) |
| 5 | In₀.₀₈Ga₀.₉₂N (400 nm) |
| 6 | In₀.₀₈Ga₀.₉₂N (30 nm, first) |
| 7 | In₀.₁₈Ga₀.₈₂N (2.6 nm, quantum well) |
| 8 | In₀.₀₈Ga₀.₉₂N (30 nm, top) |
| 100-199 | Interface layers (ID = 100 + interface number) |
| 999 | Air (if `etch_mode = 'air_material'`) |

## Usage Example with FEM

```matlab
% Generate mesh
opts.target_dz = 15;
opts.sfName = 'H27';
mesh = generatePillarQuarterMesh(opts);

% Setup finite element analysis
sf = mesh.sf;  % Use the shape function from mesh
fe = SolidElasticElem(sf, mesh.elems);

% Assign materials based on mesh.matID
for i = 1:max(mesh.matID)
    elemIndices = find(mesh.matID == i);
    % Create and assign material for each ID
    % mat = SolidMaterial(...);
    % fe.setMaterial(mat, elemIndices);
end

% Continue with analysis setup...
```

## Output Verification

When you run `generatePillarQuarterMesh()`, it prints:

1. **Configuration**: Parameters used
2. **Layer stack**: Number of physical and interface layers
3. **Verification checks**:
   - All interface layers have `nz = 1`
   - z=0 and zBank are exact slab boundaries
   - Maximum z outside pillar ≈ zBank
   - All node z-coordinates match zLevels (no warping)

## Testing

Run the test script to verify installation:

```matlab
cd examples/Pillar
testPillarMesh
```

This will:
- Generate meshes with different configurations
- Create visualization plots
- Verify mesh quality
- Print summary statistics

## Coordinate Reference

The mesh uses a right-handed coordinate system:
- **x-axis**: Points right (0 → Rout)
- **y-axis**: Points up in plane (0 → Rout)
- **z-axis**: Points up vertically (substrate at bottom, pillar extends upward)

The quarter model represents one quadrant of the full axisymmetric structure. To visualize the full structure, you would mirror across x=0 and y=0 planes.

## Advanced Usage

### Custom Layer Stack

To modify the layer stack, edit the `buildLayerStack()` function in `generatePillarQuarterMesh.m`. Ensure:
1. Porous layer ends exactly at z=0
2. No layer straddles z=0 or zBank
3. Interfaces are added between all adjacent layers

### Circle-to-Square Transition

The outer zone uses a smooth `smoothstep` function to warp from circular inner boundary to square outer boundary:

```matlab
blend = smoothstep(s);  % s ∈ [0,1]
x = (1-blend) * x_circle + blend * x_square
```

This ensures compatible meshing with periodic tiling.

## Known Limitations

1. The trench profile is qualitative and uses smooth analytic functions (not measured from fabrication)
2. Very fine meshes (target_dz < 5 nm) may have high node/element counts
3. The mesh assumes material is removed by etching (void), not deformed

## References

- Shape functions: See `math/ShapeFunctionH27.m` and `math/ShapeFunctionH8.m`
- Mesh utilities: See `mesh/Mesh.m`
- Example FEM models: See `examples/models/` directory

## Author

Generated for the MorphoGenVibrations framework (2026).
