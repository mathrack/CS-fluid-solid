for i in hexa tet
do
  for j in 05 10 15 20 30 40 60 80
  do
    echo $i$j
    ~/gmsh/gmsh-2.12.0-Linux/bin/gmsh -1 -2 -3 $i$j.geo > log.$i$j.geo
    ~/build-saturne-git/libexec/code_saturne/cs_preprocess -o $i$j $i$j.msh > log.$i$j.msh
  done
done
