# Robotic manipulator arm

Multi-load topology optimisation of a robotic arm, including the linked-density
formulation that ties repeated modules to one design field.

| | |
|---|---|
| **Branch of record** | `CAS_Arm` (tagged `paper/cas-arm`) |
| **Draft** | [Overleaf](https://www.overleaf.com/project/688de221e7936ceed901d5a6) — in preparation |
| **Publication** | *TODO — add the DOI once published; the Overleaf link is access-controlled and is not a citation* |

## Layout

- `ManipulatorArm_MainScript.m`, `ManipulatorArmBig_MainScript.m` — entry points
- `test*_*.m` — the individual experiments (linked/unlinked SIMP, stress intensity, buckling)
- `runReferenceModuleMultiLoadSIMP.m`, `sweepReferenceModuleSIMP.m` — reference-module studies
- `HelperFunctions/` — everything this study needs that the library should not carry
- `Results/` — captured runs

## Why this folder exists

This study is the worked example of the convention in §11 of
`docs/branch-integration-plan.md`: it moved out of the shared `examples/models/`
and `examples/topologyOpt/` directories in January 2026 and grew to 150 files.
The effect on merge cost is the argument for the convention — eighteen months of
parallel work against the `Vibrations` branch produced **one** conflicting
example script, against ten for the studies that stayed in shared directories.

`HelperFunctions/` holds five helpers demoted from `design/` during the Phase 3
triage. They are study code: every caller is in this folder. The rule is that a
capability earns its place in the library when a *second* study wants it, not
before.

## Library capabilities this study depends on

`StressIntensityMulti{,Av,Max}TopologyOptimization`, `Frame3D`, and the
multi-RHS load API (`createNextRightHandSideVector`,
`setCurrentLoadToRightHandSideVectors`). See decision D3 in
`docs/integration-rules.md`: this branch refactored the `Multi*` family behind a
strategy-parameterised base, and the merge keeps that structure while restoring
the per-load-case weighting the refactor dropped.
