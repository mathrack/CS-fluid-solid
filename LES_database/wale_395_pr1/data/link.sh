set -e
fid='00001'
scalar[1]='diric'
scalar[2]='neuma'
scalar[3]='11'
scalar[4]='12'
scalar[5]='13'
file='../raw/results'
ln -s ${file}.u_mean.${fid}.x.1d.dat um1d.dat
ln -s ${file}.u_mean.${fid}.y.1d.dat vm1d.dat
ln -s ${file}.u_mean.${fid}.z.1d.dat wm1d.dat
ln -s ${file}.rij.${fid}.rxx.1d.dat uum1d.dat
ln -s ${file}.rij.${fid}.ryy.1d.dat vvm1d.dat
ln -s ${file}.rij.${fid}.rzz.1d.dat wwm1d.dat
ln -s ${file}.rij.${fid}.rxy.1d.dat uvm1d.dat
ln -s ${file}.rij.${fid}.rxz.1d.dat uwm1d.dat
ln -s ${file}.rij.${fid}.ryz.1d.dat vwm1d.dat

ln -s ${file}.u_grad.${fid}.r11.1d.dat dudxm1d.dat
ln -s ${file}.u_grad.${fid}.r12.1d.dat dudym1d.dat
ln -s ${file}.u_grad.${fid}.r13.1d.dat dudzm1d.dat
ln -s ${file}.u_grad.${fid}.r21.1d.dat dvdxm1d.dat
ln -s ${file}.u_grad.${fid}.r22.1d.dat dvdym1d.dat
ln -s ${file}.u_grad.${fid}.r23.1d.dat dvdzm1d.dat
ln -s ${file}.u_grad.${fid}.r31.1d.dat dwdxm1d.dat
ln -s ${file}.u_grad.${fid}.r32.1d.dat dwdym1d.dat
ln -s ${file}.u_grad.${fid}.r33.1d.dat dwdzm1d.dat

ln -s ${file}.mu_t.${fid}.1d.dat mutm1d.dat

ln -s ${file}.u_diss.${fid}.r11.1d.dat dudx2m1d.dat
ln -s ${file}.u_diss.${fid}.r12.1d.dat dudy2m1d.dat
ln -s ${file}.u_diss.${fid}.r13.1d.dat dudz2m1d.dat
ln -s ${file}.u_diss.${fid}.r21.1d.dat dvdx2m1d.dat
ln -s ${file}.u_diss.${fid}.r22.1d.dat dvdy2m1d.dat
ln -s ${file}.u_diss.${fid}.r23.1d.dat dvdz2m1d.dat
ln -s ${file}.u_diss.${fid}.r31.1d.dat dwdx2m1d.dat
ln -s ${file}.u_diss.${fid}.r32.1d.dat dwdy2m1d.dat
ln -s ${file}.u_diss.${fid}.r33.1d.dat dwdz2m1d.dat

ln -s ${file}.mu_t_u_grad.${fid}.r11.1d.dat mududxm1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r12.1d.dat mududym1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r13.1d.dat mududzm1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r21.1d.dat mudvdxm1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r22.1d.dat mudvdym1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r23.1d.dat mudvdzm1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r31.1d.dat mudwdxm1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r32.1d.dat mudwdym1d.dat
ln -s ${file}.mu_t_u_grad.${fid}.r33.1d.dat mudwdzm1d.dat

ln -s ${file}.mu_t_u_diss.${fid}.1d.dat mutudissm1d.dat

# Scalars
for i in "${scalar[@]}"
do
  ln -s ${file}.${i}_mean.${fid}.1d.dat ${i}_mean.dat
  ln -s ${file}.${i}_variance.${fid}.1d.dat ${i}_var.dat

  ln -s ${file}.${i}_flux.${fid}.x.1d.dat ${i}_ux.dat
  ln -s ${file}.${i}_flux.${fid}.y.1d.dat ${i}_uy.dat
  ln -s ${file}.${i}_flux.${fid}.z.1d.dat ${i}_uz.dat

  ln -s ${file}.${i}_grad.${fid}.x.1d.dat ${i}_dtdx.dat
  ln -s ${file}.${i}_grad.${fid}.y.1d.dat ${i}_dtdy.dat
  ln -s ${file}.${i}_grad.${fid}.z.1d.dat ${i}_dtdz.dat

  ln -s ${file}.${i}_diss.${fid}.rxx.1d.dat ${i}_dtdx2.dat
  ln -s ${file}.${i}_diss.${fid}.ryy.1d.dat ${i}_dtdy2.dat
  ln -s ${file}.${i}_diss.${fid}.rzz.1d.dat ${i}_dtdz2.dat
  ln -s ${file}.${i}_diss.${fid}.rxy.1d.dat ${i}_dtdxy.dat
  ln -s ${file}.${i}_diss.${fid}.rxz.1d.dat ${i}_dtdxz.dat
  ln -s ${file}.${i}_diss.${fid}.ryz.1d.dat ${i}_dtdyz.dat

  ln -s ${file}.${i}_mu_t_t_grad.${fid}.x.1d.dat ${i}_mu_dtdx.dat
  ln -s ${file}.${i}_mu_t_t_grad.${fid}.y.1d.dat ${i}_mu_dtdy.dat
  ln -s ${file}.${i}_mu_t_t_grad.${fid}.z.1d.dat ${i}_mu_dtdz.dat

  ln -s ${file}.${i}_mu_t_t_diss.${fid}.1d.dat ${i}_mu_diss.dat
done
