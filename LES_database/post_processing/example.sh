#!/bin/bash -
set -e
#
# File id
#
fid='00006'
#
# Scalar names
#
scalar[1]='diric'
scalar[2]='neuma'
scalar[3]='g05a05'
scalar[4]='g05a1'
scalar[5]='g05a2'
scalar[6]='g1a05'
scalar[7]='g1a1'
scalar[8]='g1a2'
scalar[9]='g2a05'
scalar[10]='g2a1'
scalar[11]='g2a2'
#
# File names
#
file='results'
geo='.geo'
ext='.3d1d'
#
# Shortcuts
#
cc3d1d=${HOME}'/bin/import_geo'
ccvect=${HOME}'/bin/scan_geo_split_vector'
ccten6=${HOME}'/bin/scan_geo_split_tensor6'
ccten9=${HOME}'/bin/scan_geo_split_tensor'
ccstat=${HOME}'/bin/compute_1d_stats'
#
# Link to .geo and data files
#
if [ -f ./${file}${geo} ]; then
  echo "Use existing .geo file."
else
  if [ -f ../${file}${geo} ]; then
    ln -s ../$file$geo
    echo "Found and link .geo file."
  else
    echo ".geo file not found. Abort."
    exit 1
  fi
fi
if [ -f ./${file}.u_mean.${fid} ]; then
  echo "Use existing data files."
else
  if [ -f ../${file}.u_mean.${fid} ]; then
    for i in ../$file.*.$fid; do ln -s $i; done
    echo "Found and link data files."
  else
    echo "Data files not found. Abort."
    exit 1
  fi
fi
#
# Generate 3d->1d mapping
#
$cc3d1d $file$geo
#
# Explode arrays in 1d elements
#
$ccten9 $file$geo $file.u_grad.* $file.u_diss.* $file.mu_t_u_grad.*
$ccten6 $file$geo $file.rij.*
$ccvect $file$geo $file.u_mean.*
for i in "${scalar[@]}"
do
  $ccten6 $file$geo $file.${i}_diss.*
  $ccvect $file$geo $file.${i}_flux.* $file.${i}_grad.* $file.${i}_mu_t_t_grad.*
done
#
# Perform spatial averaging
#
$ccstat $file$geo$ext $file.mu_t.* $file.mu_t_u_diss.*
for i in "${scalar[@]}"
do
  $ccstat $file$geo$ext $file.${i}_mean.* $file.${i}_variance.* $file.${i}_mu_t_t_diss.*
done
$ccstat $file$geo$ext $file.*.$fid.*[1-3x-z]
