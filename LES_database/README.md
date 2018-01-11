The present folder contains the data and source code used to perform LES of the turbulent channel flow with conjugate heat transfer using *Code_Saturne*.

The folders with a name starting with *mesh_* contain all the scripts and source code to generate the fluid and solid meshes of the turbulent channel flow for the various Reynolds number investigated here.

The folder *post_processing* contains some source code used to post-process our LES. The combination of the resulting programs allows us to read the binary Ensight files produced by Code_Saturne, average it in the homogeneous directions and output the corresponding 1D profile in a text file.

The folders with a name starting with *wale_* contain all the source code and results for our LES of the turbulent channel flow with conjugate heat transfer at the given Reynolds and Prandtl numbers.
