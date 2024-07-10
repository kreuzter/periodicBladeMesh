# periodicBladeMesh
(py)gmsh based tool for the simpliest turbomachinery mesh

## Test cases:

### [SPLEEN](https://doi.org/10.5281/zenodo.7264761) turbine blade

- LE and TE definition included in the suction and pressure side
- two setups available:
  - extruded BL (gmsh internal tool)
  - transfinite mesh around blade (theoretically more stable and versatile)
- periodicities matching in geometry and mesh (the solver does not need to perform any kind of interpolation)

### [schreiber84](https://doi.org/10.1115/1.3239561) transonic compressor blade

- LE and TE computed from chord, suction and pressure side with assumption of elliptical edge
- extruded BL
- periodicities matching in geometry
