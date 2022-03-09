# SSMTool 2.2: Computation of invariant manifolds in high-dimensional mechanics problems
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4614201.svg)](https://doi.org/10.5281/zenodo.4614201)

This package computes invariant manifolds in high-dimensional dynamical systems using the *Parametrization Method* with special attention to the computation of Spectral Submanifolds (SSM) and forced response curves in finite element models.

These invariant manifolds are computed in the physical coordinates using only the master modes resulting in efficient and feasible computations for high-dimensional finite-element problems. Additionally, the user has an option to choose among the graph or normal form style of parametrization. The computational methodology is described in the following article:

[1] **Jain, S. & Haller, G. How to compute invariant manifolds and their reduced dynamics in high-dimensional finite element models. Nonlinear Dyn (2021). https://doi.org/10.1007/s11071-021-06957-4**

The theoretical and computational aspects for analyzing systems with internal resonances via multi-dimensional SSMs are given in the following articles:

[2] **Li, M., Jain, S.  & Haller, G. Nonlinear analysis of forced mechanical systems with internal resonance using spectral submanifolds -- Part I: Periodic response and forced response curve. Preprint (2021) https://arxiv.org/abs/2106.05162**

[3] **Li, M & Haller, G. Nonlinear analysis of forced mechanical systems with internal resonance using spectral submanifolds -- Part II: Bifurcation and quasi-periodic response. Preprint (2021) https://arxiv.org/abs/2108.08152**

In this version, we demonstrate the computational methodology over the following small academic examples as well high-dimensional finite element problems using the FE package [4]

First-order examples: 
- BenchmarkSSM1stOrder: computation of 1D stable SSM of a two-dimensional system
- Lorenz1stOrder: computation of 1D unstable SSM of the Lorenz system
- CharneyDeVore1stOrder: computation of 1D/2D unstable SSMs of a six-dimensional system
- Complex Dyn: example of 4D first-order dynamical system with complex coefficients benchmarked against SSMTool 1.0.

Second-order examples:
- Oscillator chain: two coupled oscillators with 1:2 internal resonance, three coupled oscillators with 1:1:1 internal resonance and n coupled oscillators without any internal resonances.
- Bernoulli beam: modeled using linear finite elements with localized nonlinearity in the form of a cubic spring with and without 1:3 **internal resonances** (IR) and demonstration of bifurcation to quasiperiodic response on a 3D torus.
- von Karman straight beam in 2D: geometrically nonlinear finite element model with and without internal resonance (1:3) and demonstration of bifurcation to quasiperiodic response on a 2D torus.
- von Karman plate in 3D: geometrically nonlinear finite element model of a square flat plate with demonstration of **parallel computing**, 1:1 internal resonance and bifurcation to quasiperiodic response on a 2D torus.
- von Karman shell-based shallow curved panel in 3D: geometrically nonlinear finite element model with and without 1:2 internal resonance.
- Prismatic beam: nonlinear beam PDE discretized using Galerkin method onto a given number of modes, comparison with the method of multiple scales, demonstration of 1:3 internal internal resonance
- AxialMovingBeam: an axially moving beam with gyroscopic and nonlinear damping forces
- TimoshenkoBeam: a cantilever Timoshenko beam carrying a lumped mass. This example demonstrates the effectiveness of SSM reduction for systems undergoing large deformations.
- NACA airfoil based aircraft wing model: shell-based nonlinear finite element model containing more than 100,000 degrees of freedom.



This package uses the following external open-source packages:

1. Continuation core (coco) https://sourceforge.net/projects/cocotools/
2. Sandia tensor toolbox: https://gitlab.com/tensors/tensor_toolbox
3. Combinator: https://www.mathworks.com/matlabcentral/fileexchange/24325-combinator-combinations-and-permutations
4. YetAnotherFECode: Zenodo http://doi.org/10.5281/zenodo.4011281
5. Tor: https://github.com/mingwu-li/torus_collocation

In order to install the program, simply run the install.m file in the main folder. The examples can be found in the examples directory.
Note: When running the examples in the livescript files (workbooks), please ensure that the MATLAB 'Current Folder' is the directory of the specific example.

Please report any issues/bugs to Shobhit Jain shjain@ethz.ch
