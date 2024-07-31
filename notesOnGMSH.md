# GMSH

gmsh is very well documented in 
[official documentation](https://gmsh.info/doc/texinfo/gmsh.html), but info useful for setup
of this tool is pointed out here.

## Algorithms

Saved in: `General.OptionsFileName`

**`Mesh.Algorithm`**: 2D mesh algorithm, options:

  - `1`: MeshAdapt, 
  - `2`: Automatic, 
  - `3`: Initial mesh only, 
  - `5`: Delaunay, 
  - `6`: Frontal-Delaunay (*default*),
  - `7`: BAMG, 
  - `8`: Frontal-Delaunay for Quads, 
  - `9`: Packing of Parallelograms, 
  - `11`: Quasi-structured Quad

**`Mesh.RecombinationAlgorithm`**: mesh recombination algorithm 
  - `0`: simple,
  - `1`: blossom (*default*),
  - `2`: simple full-quad,
  - `3`: blossom full-quad

## Available formats

- `1`: msh,
- `2`: unv, 
- `10`: auto, 
- `16`: vtk, 
- `19`: vrml, 
- `21`: mail, 
- `26`: pos stat, 
- `27`: stl, 
- `28`: p3d, 
- `30`: mesh, 
- `31`: bdf, 
- `32`: cgns, 
- `33`: med, 
- `34`: diff, 
- `38`: ir3, 
- `39`: inp, 
- `40`: ply2, 
- `41`: celum, 
- `42`: su2, 
- `47`: tochnog, 
- `49`: neu, 
- `50`: matlab