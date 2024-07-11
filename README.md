# periodicBladeMesh
(py)gmsh based tool for the simpliest turbomachinery mesh

Internal references can be used in json files. E.g.:
```
  "scaled by"  : {"$ref": "#/geometry/true chord"},
```

## Test cases:

### [SPLEEN](https://doi.org/10.5281/zenodo.7264761) turbine blade

- LE and TE definition included in the suction and pressure side
- two setups available:
  - extruded BL (gmsh internal tool)
  - transfinite mesh around blade (theoretically more stable and versatile)
- periodicities matching in geometry and mesh (the solver does not need to perform any kind of interpolation)
