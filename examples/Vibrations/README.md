# Beam vibrations

Free-vibration and self-weight topology optimisation of beams and a building
frame, plus the buckling work this study grew out of.

| | |
|---|---|
| **Branch of record** | `Vibrations` (tagged `paper/vibrations`) |
| **Buckling ancestor** | `Buckling` (tagged `paper/buckling`) — zero unique commits, folded in here |
| **Publication** | *TODO — add the DOI once published* |

## Layout

- `BeamFromPaper*.m`, `BeamClampedFromPaper.m` — the beam cases reproduced from the paper
- `BeamSelfWeightTopOpt*.m`, `BuildingSelfWeightTopOpt.m` — self-weight topology optimisation
- `Buckling/` — the buckling study: 2D and 3D cantilevers, and the manipulator case
- `*.json` — run configurations

## Library capabilities this study depends on

`LinearNaturalVibration`, `LinearStability`, `ElasticHarmonicVibrations`,
`SecondOrderElasticityWeighted`, and per-load-case weighting (`alphas`) in the
`StressIntensityMulti*` optimisers — see decision D3 in
`docs/integration-rules.md`, which exists because a naive merge would have
dropped that weighting and left the sweep running with every weight equal.
