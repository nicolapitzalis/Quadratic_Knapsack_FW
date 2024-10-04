# Frank-Wolfe Algorithm Variants Comparison
This work was completed as part of the Computational Mathematics for Learning and Data Analysis course at the University of Pisa.
This project aims to compare the performance of three variants of the Frank-Wolfe (FW) algorithm—`fw_og`, `fw_away_step`, and `fw_pairwise`—in solving convex quadratic programming problems. The project runs a series of experiments across different configurations and evaluates the computational time and convergence behavior of each algorithm.

## Overview

The project includes the following components:
1. **Algorithms**:
   - `fw_og`: The original Frank-Wolfe algorithm.
   - `fw_away_step`: The Away-Step variant of Frank-Wolfe.
   - `fw_pairwise`: The Pairwise variant of Frank-Wolfe.

2. **Experiments**:
   - Conducts experiments by varying problem parameters such as:
     - Problem size (`n`)
     - Constraint activity (`actv`)
     - Matrix density (`density`)
     - Reciprocal condition number (`rc`)
   - Compares the performance of the algorithms in terms of computational time and convergence (duality gap).

3. **Visualization**:
   - Generates comparison plots of time and convergence for each experiment.
   - The plots show how the different algorithm variants perform under various conditions, such as varying problem size or condition number.

## Structure

- **Main Script**:
  - `run_experiments.m`: Executes the experiments by varying problem parameters and saves the results.

- **Frank-Wolfe Variants**:
  - `fw_og.m`: Original Frank-Wolfe algorithm.
  - `fw_away_step.m`: Away-Step variant of Frank-Wolfe.
  - `fw_pairwise.m`: Pairwise variant of Frank-Wolfe.

- **Helper Functions**:
  - `test_quad.m`: Runs each Frank-Wolfe variant for a given configuration.
  - `check_feasibility.m`: Verifies the feasibility of solutions with respect to the constraints.
  - `compute_pyramidal_width.m`: Approximates the pyramidal width of the polytope.
  - `curvature_function.m`: Computes the curvature of the objective function.
  - `get_gap.m`: Calculates the duality gap for a given solution.
  - `compute_diameter_polytope.m`: Computes the diameter of the polytope defined by the problem's bounds.

## Run Experiments

1. **Run Experiments**:
   To run the predefined experiments, execute:
   ```matlab
   run_experiments();
   ```

2. **View Results**:
   After running the experiments, view the comparison plots saved in the `results/` directory.

3. **Customize Parameters**:
   You can modify the `run_experiments.m` file to change the problem parameters such as `n`, `actv`, `density`, or `rc` to run different experiments.

## Test Frank-Wolfe Variants

To test the Frank-Wolfe variants with a specific configuration, you can use the `test_quad` function. For example:
```matlab
test_quad(alg, epsilon, max_iter, n, actv, density, rc, trials, bound_flag, path)
```
where:
- `alg`: The optimization algorithm to test, one of: 'fw_og', 'fw_away_step', 'fw_pairwise', or 'quadprog'.
- `epsilon`: The convergence tolerance of the duality gap for the Frank-Wolfe variants.
- `max_iter`: Maximum number of iterations allowed for the algorithms.
- `n`: Problem size.
- `actv`: Fraction of active constraints in the problem.
- `density`: Matrix density.
- `rc`: Reciprocal condition number.
- `trials`: Number of trials to run.
- `bound_flag`: If true, plots the dual gap along with upper bounds during the experiment.
- `path`: Path to save the results.

## Requirements

- MATLAB with Optimization Toolbox (for `linprog`, `fmincon`, etc.)
- The experiments are set up to run in parallel if available (`parpool`).

## Output

- **Results**: The logs and plots for each experiment are saved in the `results/` directory, structured by the experiment type.
- **Plots**: PNG files of comparison plots for each experiment.
