# SSMTool 2.5: Computation of invariant manifolds in high-dimensional mechanics problems
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4614201.svg)](https://doi.org/10.5281/zenodo.4614201)

What's new in SSMTool 2.5
- Computation of forced response surface (FRS) via SSM-based model reduction,
- Extraction of ridges and trenches of FRS without computing the surface,
- Explicit adjoints and gradients for the optimization of periodic orbits via SSM-based ROMs. 

This package computes invariant manifolds in high-dimensional dynamical systems using the *Parametrization Method* with special attention to the computation of Spectral Submanifolds (SSM) and forced response curves in finite element models.

These invariant manifolds are computed in the physical coordinates using only the master modes resulting in efficient and feasible computations for high-dimensional finite-element problems. Additionally, the user has an option to choose among the graph or normal form style of parametrization. The computational methodology is described in the following article:

[1] **Jain, S. & Haller, G. How to compute invariant manifolds and their reduced dynamics in high-dimensional finite element models. Nonlinear Dyn (2021). https://doi.org/10.1007/s11071-021-06957-4**

The theoretical and computational aspects for analyzing systems with internal resonances via multi-dimensional SSMs are given in the following articles:

[2] **Li, M., Jain, S.  & Haller, G. Nonlinear analysis of forced mechanical systems with internal resonance using spectral submanifolds, Part I: Periodic response and forced response curve. Nonlinear Dyn 110, 1005–1043 (2022). https://doi.org/10.1007/s11071-022-07714-x**

[3] **Li, M & Haller, G. Nonlinear analysis of forced mechanical systems with internal resonance using spectral submanifolds, Part II: Bifurcation and quasi-periodic response. Nonlinear Dyn 110, 1045–1080 (2022). https://doi.org/10.1007/s11071-022-07476-6**

How SSMs are extended to constrained mechanical systems are discussed in the following article:

[4] **Li, M., Jain, S.  & Haller, G. Model reduction for constrained mechanical systems via spectral submanifolds. Nonlinear Dyn 111(10): 8881-8911 (2023).
https://doi.org/10.1007/s11071-023-08300-5**

The treatment of systems subject to parametric resonance via higher-order approximations of nonautonomous SSMs is described in the following article:

[5] **Thurnher, T., Haller, G.  & Jain, S. Nonautonomous spectral submanifolds for model reduction of nonlinear mechanical systems under parametric resonances. Preprint (2023). Available on arXiv: 
https://doi.org/10.48550/arXiv.2307.10240**

The use of SSM-based ROMs to extract forced response surfaces (FRSs) and their ridges and trenches via parameter continuation is discussed in the following article:

[6] **Li, M., Jain, S.  & Haller, G. Fast computation and characterization of forced response surface via spectral submanifolds and parameter continuation. Preprint (2023). Available on arXiv: 
https://doi.org/10.48550/arXiv.2310.14850**

In this version, we demonstrate the computational methodology over the following small academic examples as well high-dimensional finite element problems using the FE package *YetAnotherFECode*

First-order examples:
- BenchmarkSSM1stOrder: computation of 1D stable SSM of a two-dimensional system
- Lorenz1stOrder: computation of 1D unstable SSM of the Lorenz system
- CharneyDeVore1stOrder: computation of 1D/2D unstable SSMs of a six-dimensional system
- Complex Dyn: example of 4D first-order dynamical system with complex coefficients benchmarked against SSMTool 1.0.

Second-order examples:
- Oscillator chain: two coupled oscillators with 1:2 internal resonance, three coupled oscillators with 1:1:1 internal resonance and n coupled oscillators without any internal resonances. [1]
- Bernoulli beam: modeled using linear finite elements with localized nonlinearity in the form of a cubic spring with and without 1:3 **internal resonances** (IR) and demonstration of bifurcation to quasiperiodic response on a 3D torus. [3]
- von Karman straight beam in 2D: geometrically nonlinear finite element model with and without internal resonance (1:3) and demonstration of bifurcation to quasiperiodic response on a 2D torus. [1,2,3]
- von Karman plate in 3D: geometrically nonlinear finite element model of a square flat plate with demonstration of **parallel computing**, 1:1 internal resonance and bifurcation to quasiperiodic response on a 2D torus. [2,3]
- von Karman shell-based shallow curved panel in 3D: geometrically nonlinear finite element model with and without 1:2 internal resonance. [2]
- Prismatic beam: nonlinear beam PDE discretized using Galerkin method onto a given number of modes, comparison with the method of multiple scales, demonstration of 1:3 internal internal resonance
- AxialMovingBeam: an axially moving beam with gyroscopic and nonlinear damping forces [2]
- PipeConveyingFluid: an **fluid-structure interaction** system with flow-induced gyroscopic and follower forces, demonstration of **asymmetric** damping and stiffness matrices and **global bifurcation**. More details can be found in https://doi.org/10.1016/j.ymssp.2022.109993
- TimoshenkoBeam: a cantilever Timoshenko beam carrying a lumped mass. This example demonstrates the effectiveness of SSM reduction for systems undergoing large deformations. [2]
- NACA airfoil based aircraft wing model: shell-based nonlinear finite element model containing more than 100,000 degrees of freedom. [1]

Constrained mechanical systems [4]
 - 3D oscillator constrained in a surface
 - Pendulums in Cartesian coordinates: single pendulum, pendulum-slider with 1:3 internal resonance, and a chain of pendulums
 - Frequency divider: two geometrically nonlinear beams connected via a revolute joint

Computation of stability diagrams and forced response curves in mechanical systems under parametric resonance: 
 - 1-DOF Mathieu oscillator with periodically-varying linear and cubic spring 
 - Coupled system of two Mathieu oscillators with cubic nonlinearity [5]
 - **Self-excited** 2-DOF-oscillator system with cubic nonlinearity, time-varying stiffness as well as external periodic forcing, where the forced response exhibits **multiple isolas**
 - Self-excited 2-DOF-oscillator system with external and parametric excitation  [5]
 - **Parametric amplification** in coupled-oscillators due to nonlinear damping and time-varying stiffness and external periodic forcing with a phase lag relative to stiffness.
 - Prismatic beam: nonlinear beam PDE discretized using Galerkin method onto a given number of modes, demonstration of axial stretching leading to parametric excitation of transverse degrees of freedom. [5]
 - Bernoulli beam attached to a spring time-varying linear stiffness along with a nonlinear damper and spring: Stability diagrams as well as forced response curves (exhibiting isolas) are constructed over two-dimensional SSMs. [5]

Forced response surface and its ridges and trenches [6]
 - Bernoulli beam with cubic nonlinear support: analytic prediction of FRS, which automatically uncovers isolas of forced response curves
 - von Karman plate with 1:1 internal resonance: extraction of FRS via two-dimensional continuation; extraction of ridges and trenches via successive parameter continuation
 - von Karman shell with 1:2 internal resonance: extraction of FRS and its ridges and trenches.

This package uses the following external open-source packages:

1. Continuation core (coco) https://sourceforge.net/projects/cocotools/
2. Sandia tensor toolbox: https://gitlab.com/tensors/tensor_toolbox
3. Combinator: https://www.mathworks.com/matlabcentral/fileexchange/24325-combinator-combinations-and-permutations
4. YetAnotherFECode: Zenodo http://doi.org/10.5281/zenodo.4011281
5. Tor: https://github.com/mingwu-li/torus_collocation

In order to install the program, simply run the install.m file in the main folder. The examples can be found in the examples directory.
Note: When running the examples in the livescript files (workbooks), please ensure that the MATLAB 'Current Folder' is the directory of the specific example.

Please report any issues/bugs to Shobhit Jain (shobhit.jain@tudelft.nl) and Mingwu Li (limw@sustech.edu.cn)
