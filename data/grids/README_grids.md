# Grid

Valentin Dall'Alba / February 2020

## Description :

The grids are note stored on git hub due to the size of the files. For access to the different grids please send an email to *valentin.dallalba-arnau@unine.ch*

The different grids and variables are the Roussillon grids, the 3D trend map and the 3D rotation maps.

## Files : 

* _roussillon_grids.gslib_ is the grid of the Continental Pliocene (PC) layer. This gslib file is composed of two variables. The first one corresponds to the PC layer created from the top and bottom topography maps of the PC. The second variable is the simulation grid which corresponds to the PC grid transformed/flattened in respect of the bottom topography of the PC layer. The transformed grid allows to simulate in 2D layer sharing similar age of deposition. Details on the grid creation process can be found in the associated jupyter notebook (*jupyter/grid_creation/createGrid.ipynb*).

* _trend_map.gslib_ corresponds to the auxiliary variable of the simulation grid. The creation process of this trend map is described in the article but details regarding the code used to simulate the variable can be found in the jupyter notebook associated (_jupyter/trend_map_creation/createTrendMap.ipynb_).

* _rotation_maps.gslib_ is the rotation maps used for rotate the pattern during the simulations. Two rotation maps were used to create a 20° of tolerance in the rotation value selected by the algorithm, details can be found in the associated jupyter notebook (*jupyter/rotation_map_creation/createRotationMap.ipynb*).



 

