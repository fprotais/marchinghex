# marchinghex
Robust hexahedral meshing using Dhondt cut approach.

Technical details are available in a technical report: [https://inria.hal.science/hal-04913435v1](https://inria.hal.science/hal-04913435v1)

The code take as an input ***a grid*** (as an hex mesh) and ***a domain***  (as a tet mesh) and extract an hexahedral mesh of the domain.
To match with the boundary, the code reads triangles and edges from the domain mesh. 

Starting from a regular grid, elements are guaranteed to have positive scaled-jacobian. 
To match the boundary, we slowly move points, making sure no elements are inverted. 
**Starting from a regular grid, we are guaranteed to obtain an hexahedral mesh free of inverted elements (with scaled-jacobian \> 0)**. 

The smoothing is done using an elliptic energy, and is currently quite slow (5000 nodes/s). I am working to improve this.

The code rely on lib ultimaille and uses its io. Supported mesh files are *geogram*, *medit* and *vtk*. 
You can display hexahedral mesh quality using [hexalab.net](https://www.hexalab.net/). 
For generating and visualising meshes in general, I recommand [gmsh](http://gmsh.info/) or [Graphite](http://alice.loria.fr/index.php?option=com_content&view=article&id=22).

**Following ultimaille, this code is under GNU Affero General Public License v3.0 (see LICENSE).**
If you are faithful to AGPL but still want to include this code in a commercial software, feel free to contact me, I am always happy to work with new people. 
And the same goes for non-commercial and academic software of course. 

# Use CMake to build the project:
```sh
git clone --recurse-submodules https://github.com/fprotais/marchinghex &&
cd marchinghex &&
mkdir build &&
cd build &&
cmake .. &&
make -j 
```
Build on windows works only with static libs. 

# Grid generation

We provided a small executable to generate input grid mesh. 
The grid generated is a regular grid with elements similar to the input tetrahedra. 
A scale change the coarsness of the grid (5 : very fine, 0.1: very coarse).

```sh
./gridgenerator ../meshes/fertility.mesh grid.mesh 1.
```
![box](https://raw.githubusercontent.com/fprotais/marchinghex/main/images/mesh_in_grid.jpg)

/!\ Note that for hexmeshing to work the best, each point of grid should ***absolutly not*** be on the boundary of the domain. 
The given grid generator is trivial, and doesn't account for this. 

# Hex meshing

The code reads on the domain both triangles and edges. 
Triangle will be use for boundary projections, as edges for features projections.
It is not guaranteed to match perfectly the boundary, as it must keep element validity, but it works pretty well.

![box](https://raw.githubusercontent.com/fprotais/marchinghex/main/images/hexmeshing.jpg)
To generate this result:
```sh
./gridgenerator ../meshes/fertility.mesh grid.mesh 1.
./marchinghex_hexmeshing grid.mesh ../meshes/fertility.mesh hexmesh.mesh
```
![box](https://raw.githubusercontent.com/fprotais/marchinghex/main/images/caohexmeshing.jpg)
To generate this result:
```sh
./gridgenerator ../meshes/CAD0.mesh grid.mesh 2.
./marchinghex_hexmeshing grid.mesh ../meshes/CAD0.mesh hexmesh.mesh
```

# Bi-material meshing

The pattern used allows for bimaterial meshing, with inside and outside mesh being compatible. 
![box](https://raw.githubusercontent.com/fprotais/marchinghex/main/images/bimaterial.jpg)
To generate this result:
```sh
./gridgenerator ../meshes/fertility.mesh grid.mesh 1.
./marchinghex_bimaterial grid.mesh ../meshes/fertility.mesh bimaterial.mesh inside.mesh outside.mesh
```

Elements produced are of good quality:
![box](https://raw.githubusercontent.com/fprotais/marchinghex/main/images/bimaterial_sj.jpg)


# Quality of hexahedra

We use a method which generates polyhedra for each cube of the grid, and then obtains hexahedra through midpoint splitting. 
"different_configurations.mesh" contains all the configuration generated with our method:
![box](https://raw.githubusercontent.com/fprotais/marchinghex/main/images/configurations.png)

You can notice that some elements have negative scaled-jacobian. 
To avoid those, \[1\] proposes a clever splitting, that gives "splitted_configurations.mesh" and "splitted_configurations_compl.mesh", each configuration and its complementary:
![box](https://raw.githubusercontent.com/fprotais/marchinghex/main/images/with_split.jpg)
To generate those results: 
```sh
./make_examples
```
*Generates the complete mesh, and subdirectories containing all differents configurations in individual meshes*






