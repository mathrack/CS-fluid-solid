The present folder contains some source code used to post-process our LES. Here, we describe only the main files. It seems important to stress that the resulting programs are not specifically robust and were not designed to be generic. They were fully tested and should work with the LES performed in this study.

* import_geo.f90: Resulting program reads the 3D *.geo* file and extract the 1D wall-normal grid and the connectivity between the 3D mesh and the 1D grid. The 1D grid and the 3D-1D mapping are stored in a specific file which has the *.3d1d* extension. The 3D binary mesh must be passed as a command-line argument.
* scan_geo_split_vector.f90, scan_geo_split_tensor6.f90, scan_geo_split_tensor.f90: Resulting programs will split vectors or tensors given as argument into 1D components. Each component is written in a new binary file. The first command-line argument must be the 3D *.geo* file. The following arguments should be the vectors or tensors to split.
* compute_1d_stats.f90: Resulting program will average scalar fields over homogeneous directions and output the 1D profile in a text file. The first command-line argument must be the *.3d1d* file obtained with import_geo. The following arguments should be the 1D scalar fields to average.

The bash script *example.sh* illustrates the combined usage of those programs.
