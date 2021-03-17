# SSMTool-2.0: Computation of invariant manifolds in high-dimensional mechanics problems

This package computes invariant manifolds in high-dimensional dynamical systems using the Parametrization Method with special attention to the computation of Spectral Submanifolds (SSM) and forced response curves in finite element models. 

Such manifolds are computed in the physical coordinates using only the master modes resulting in efficient and feasible computations for high-dimensional problems. Additionally, the user has an option to choose among the graph or normal form style of parametrization. 

In this pre-release version, we demonstrate the computational methodology over small academic examples as well high-dimensional finite element problems using the FE package [4].

In this version, we have included a demonstration of SSM computation over the following finite element examples

- von Karman straight beam in 2D
- von Karman shell-based shallow curved panel in 3D
- NACA airfoil based aircraft wing model

This package uses the following external open-source packages:

1. Continuation core (coco) https://sourceforge.net/projects/cocotools/
2. Sandia tensor toolbox: https://gitlab.com/tensors/tensor_toolbox
3. Combinator: https://www.mathworks.com/matlabcentral/fileexchange/24325-combinator-combinations-and-permutations
4. YetAnotherFECode: Zenodo http://doi.org/10.5281/zenodo.4011282

In order to install the program, simply run the install.m file in the main folder. The examples can be found in the examples directory.
Note: When running the examples in the livescript files (workbooks), please ensure that the MATLAB 'Current Folder' is the directory of the specific example.
