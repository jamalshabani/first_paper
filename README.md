# Numerical Experiments for Compliant Morphing Structures

This repository contains Jamal Shabani's numerical code used to produce results for the paper:

**Systematic Design of Compliant Morphing Structures: A Phase-Field Approach**

It also contains computations associated with the doctoral thesis:

**Systematic design of compliant morphing structures with stimulus as design and state variable**

The implementation uses Firedrake for finite-element computations and PETSc/TAO for numerical optimization.

## Research context

The research concerns the systematic design of compliant structures containing structural and responsive materials. The numerical experiments investigate material distributions and stimulus-driven deformation using a phase-field approach.

This repository focuses on the computations associated with the paper and thesis. Broader FireTop experiments and ongoing trajectory work are maintained separately.

## Related repositories

| Repository | Purpose |
| --- | --- |
| [optimal](https://github.com/jamalshabani/optimal) | Main FireTop research code. |
| [first_paper](https://github.com/jamalshabani/first_paper) | Paper and thesis computations. |
| [trajectory](https://github.com/jamalshabani/trajectory) | Ongoing trajectory topology optimization research. |

## Experiment layout

| Source | Experiment family | Default mesh |
| --- | --- | --- |
| `directFirstComputation.py` | Direct first computation | `motion.msh` |
| `alternateFirstComputation.py` | Alternate first computation | `motion.msh` |
| `secondComputation.py` | Second computation | `motion.msh` |
| `hexagonal.py` | Hexagonal-domain computation | `hexagonal.msh` |

The repository also contains:

- `RunScripts/`: additional experiment-launching material.
- `directFirstComputationRatio*`: archived direct-computation results.
- `alternateFirstComputationRatio*`: archived alternate-computation results.
- `secondComputationRatio*`: archived second-computation results.
- `hexagonalComputationRatio*`: archived hexagonal-domain results.

Use the actual parameter values and commands to identify a material contrast. Directory names alone do not fully specify an experiment.

## Dependencies

The numerical solvers require:

- Python 3.
- Firedrake.
- PETSc with TAO, accessed through `petsc4py`.
- NumPy.

The paper scripts explicitly import `VTKFile` from `firedrake.output`.

Check the required imports inside your Firedrake environment:

```bash
python3 -c "import firedrake, numpy; from firedrake.output import VTKFile; from petsc4py import PETSc; print(PETSc.Sys.getVersion())"
```

A fully pinned historical environment is not specified here. Record dependency versions when reproducing results.

## Running an example

Run commands from the repository root.

Inspect the solver options:

```bash
python3 directFirstComputation.py --help
```

Run an initial example:

```bash
mkdir -p runs/direct-example

python3 directFirstComputation.py \
  -m motion.msh \
  -o runs/direct-example/beam.pvd \
  -tao_monitor \
  -tao_max_it 20 \
  -er 1.0 \
  -es 0.01
```

The paper solvers pass `-o` directly to `VTKFile`. Supply a **PVD filename**, including its parent directory.

The reduced iteration limit is suitable for an initial execution check. Reproducing a converged research result requires the corresponding experiment settings and convergence criteria.

## Main parameters

| Options | Meaning |
| --- | --- |
| `-m`, `-o` | Input mesh and output PVD filename. |
| `-er`, `-es` | Responsive and structural elastic moduli. |
| `-vr`, `-vs` | Responsive and structural volume parameters. |
| `-lr`, `-ls` | Material Lagrange multiplier parameters. |
| `-k`, `-e` | Regularization weight and phase-field parameter. |
| `-p` | Interpolation exponent. |
| `-s` | Initial stimulus parameter. |
| `-tao_type` | Optimization algorithm. |
| `-tao_ls_type` | Line-search selection. |
| `-tao_max_it` | Maximum optimization iterations. |
| `-tao_max_funcs` | Maximum function evaluations. |
| `-tao_monitor`, `-tao_view` | Optimization diagnostics. |

Defaults vary between formulations. Boundary conditions, initial fields, and objectives are defined in the corresponding Python source.

## Reproducing the research results

1. Identify the experiment in the paper or thesis.
2. Select the corresponding solver and mesh.
3. Inspect the available launch scripts and parameter values.
4. Match the elastic moduli, regularization, material volume settings, stimulus, and optimizer tolerances.
5. Write results to a fresh output path.
6. Save the command, source commit, dependency versions, and solver log.
7. Compare the converged fields and reported quantities with the research result.

Retain the supplied mesh boundary labels. Changes to the mesh require reviewing the boundary conditions in the solver.

This README provides an entry point to the computations; it does not establish a complete figure-by-figure mapping. Archived outputs should be considered alongside the original input settings and convergence history.

## Viewing results

Open the generated `.pvd` collections in ParaView or another VTK-compatible viewer.

Keep the associated data files in their expected locations. Moving a PVD file without its referenced files can prevent the results from loading.

The repository contains large simulation archives, so downloading the complete repository may require substantial storage.

## Publications and citation

### Associated paper

J. Shabani, K. Bhattacharya, and B. Bourdin.

*Systematic Design of Compliant Morphing Structures: A Phase-Field Approach.*

[Read the preprint](https://arxiv.org/abs/2411.06289)

### Doctoral thesis

Jamal Shabani.

*Systematic design of compliant morphing structures with stimulus as design and state variable.*

McMaster University.

[Thesis Publication Link](https://macsphere.mcmaster.ca/items/8108f6e2-92a0-4529-a892-baa3c2526f2d)

When using these computations, cite the associated research and record the repository URL and commit used.

## Author

**Jamal Shabani**

[GitHub](https://github.com/jamalshabani)
