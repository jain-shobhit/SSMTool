GUIDELINES
- ContinuumElement.m should contain only generic methods
- specific methods (e.g. stiffness_derivative.m, which is
  used only to compute Modal Derivatives in some projection
  approaches) should be included in the @ContinuumElement
  folder as separate functions to keep the code modular and
  clean.


List of available elements:
- Quad4 (4-noded quadrilateral, 2D, linear shape functions)
- Quad8 (8-noded quadrilateral, 2D, quadratic shape functions)
- Hex8  ( 8-noded hexahedron, linear shape functions)
- Hex20 (20-noded hexahedron, quadratic shape functions)
- Tet4  ( 4-noded tetrahedron, linear shape functions)
- Tet10 (10-noded tetrahedron, quadratic shape functions)
- Wed15 (15-noded wedge, quadratic shape functions)


Last modified: 9 Dec 2021
Jacopo Marconi, PhD, Politecnico di Milano