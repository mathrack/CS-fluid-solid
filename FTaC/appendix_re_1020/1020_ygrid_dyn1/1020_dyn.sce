nb_point=172;
nx=200;
nz=200;
lx=6.283;
ly=2.;
lz=3.142;
mu=0.00004826;
re=1.0/mu;
dx = lx / nx;
dz = lz / nz;

maxfiles(100);
format('e',16);

y=zeros(nb_point,1);
ume =zeros(y);
uume =zeros(y);
vvme =zeros(y);
wwme =zeros(y);
uvme =zeros(y);
dudy =zeros(y);
mu_t =zeros(y);

function [yy,var] = read_my_stat_1d(filename,nb)
  [fd_um, ier]=mopen(filename,'r');
  for i=1:nb
    tmp_m=mfscanf(fd_um,'%f %f'); 
    yy(i)=1.0+tmp_m(1);
    var(i)=tmp_m(2);
  end
  mclose(fd_um);
endfunction;

[y,ume]      = read_my_stat_1d('./data/um1d.dat',nb_point);
[y,vme]      = read_my_stat_1d('./data/vm1d.dat',nb_point);
[y,wme]      = read_my_stat_1d('./data/wm1d.dat',nb_point);
[y,uume]     = read_my_stat_1d('./data/uum1d.dat',nb_point);
[y,vvme]     = read_my_stat_1d('./data/vvm1d.dat',nb_point);
[y,wwme]     = read_my_stat_1d('./data/wwm1d.dat',nb_point);
[y,uvme]     = read_my_stat_1d('./data/uvm1d.dat',nb_point);
[y,dudy]     = read_my_stat_1d('./data/dudym1d.dat',nb_point);
[y,mu_t]     = read_my_stat_1d('./data/mutm1d.dat',nb_point);

[y,ume]      = read_my_stat_1d('./data/um1d.dat',nb_point);
write_csv([y,ume,vme,wme], "./csv/raw_1.csv"," ");
write_csv([y,uume,vvme,wwme,uvme], "./csv/raw_2.csv"," ");

ny=nb_point/2;
ym=zeros(ny,1);
yp=zeros(ny,1);
um=zeros(ny,1);
up=zeros(ny,1);
mu_tm=zeros(ny,1);
mu_tp=zeros(ny,1);
uum=zeros(ny,1);
uup=zeros(ny,1);
vvm=zeros(ny,1);
vvp=zeros(ny,1);
wwm=zeros(ny,1);
wwp=zeros(ny,1);
uvm=zeros(ny,1);
uvp=zeros(ny,1);

for i=1:ny
   um(i)= ume(i);
mu_tm(i)= mu_t(i);
  uum(i)= uume(i);
  vvm(i)= vvme(i);
  wwm(i)= wwme(i);
  uvm(i)= uvme(i);

   up(i)= ume(nb_point-i+1);
mu_tp(i)= mu_t(nb_point-i+1);
  uup(i)= uume(nb_point-i+1);
  vvp(i)= vvme(nb_point-i+1);
  wwp(i)= wwme(nb_point-i+1);
  uvp(i)= uvme(nb_point-i+1);
end

utau_m=sqrt( (abs(dudy(1)))/re );
utau_p=sqrt( (abs(dudy(nb_point)))/re );
utau=sqrt(0.5*(utau_m*utau_m+utau_p*utau_p))

for i=1:ny
  ym(i)=re*utau_m*y(i);
  yp(i)=re*utau_p*(2-y(nb_point-i+1));

  um(i)=um(i)/utau_m;
  up(i)=up(i)/utau_p;

  uum(i)=uum(i)/(utau_m*utau_m);
  uup(i)=uup(i)/(utau_p*utau_p);
  vvm(i)=vvm(i)/(utau_m*utau_m);
  vvp(i)=vvp(i)/(utau_p*utau_p);
  wwm(i)=wwm(i)/(utau_m*utau_m);
  wwp(i)=wwp(i)/(utau_p*utau_p);

  uvm(i)=uvm(i)/(utau_m*utau_m);
  uvp(i)=uvp(i)/(utau_p*utau_p);
end

[y,ume]      = read_my_stat_1d('./data/um1d.dat',nb_point);
write_csv([ym,um,mu_tm/mu], "./csv/moy1.csv"," ");
write_csv([yp,up,mu_tp/mu], "./csv/moy2.csv"," ");
write_csv([ym,uum,vvm,wwm,-uvm], "./csv/fluct1.csv"," ");
write_csv([yp,uup,vvp,wwp,uvp], "./csv/fluct2.csv"," ");

re*utau

[y,dudx] = read_my_stat_1d('./data/dudxm1d.dat',nb_point);
//[y,dudy] = read_my_stat_1d('./data/dudym1d.dat',nb_point);
[y,dudz] = read_my_stat_1d('./data/dudzm1d.dat',nb_point);

[y,dvdx] = read_my_stat_1d('./data/dvdxm1d.dat',nb_point);
[y,dvdy] = read_my_stat_1d('./data/dvdym1d.dat',nb_point);
[y,dvdz] = read_my_stat_1d('./data/dvdzm1d.dat',nb_point);

[y,dwdx] = read_my_stat_1d('./data/dwdxm1d.dat',nb_point);
[y,dwdy] = read_my_stat_1d('./data/dwdym1d.dat',nb_point);
[y,dwdz] = read_my_stat_1d('./data/dwdzm1d.dat',nb_point);

[y,mududx] = read_my_stat_1d('./data/mududxm1d.dat',nb_point);
[y,mududy] = read_my_stat_1d('./data/mududym1d.dat',nb_point);
[y,mududz] = read_my_stat_1d('./data/mududzm1d.dat',nb_point);

[y,mudvdx] = read_my_stat_1d('./data/mudvdxm1d.dat',nb_point);
[y,mudvdy] = read_my_stat_1d('./data/mudvdym1d.dat',nb_point);
[y,mudvdz] = read_my_stat_1d('./data/mudvdzm1d.dat',nb_point);

[y,mudwdx] = read_my_stat_1d('./data/mudwdxm1d.dat',nb_point);
[y,mudwdy] = read_my_stat_1d('./data/mudwdym1d.dat',nb_point);
[y,mudwdz] = read_my_stat_1d('./data/mudwdzm1d.dat',nb_point);

[y,dudx2] = read_my_stat_1d('./data/dudx2m1d.dat',nb_point);
[y,dudy2] = read_my_stat_1d('./data/dudy2m1d.dat',nb_point);
[y,dudz2] = read_my_stat_1d('./data/dudz2m1d.dat',nb_point);

[y,dvdx2] = read_my_stat_1d('./data/dvdx2m1d.dat',nb_point);
[y,dvdy2] = read_my_stat_1d('./data/dvdy2m1d.dat',nb_point);
[y,dvdz2] = read_my_stat_1d('./data/dvdz2m1d.dat',nb_point);

[y,dwdx2] = read_my_stat_1d('./data/dwdx2m1d.dat',nb_point);
[y,dwdy2] = read_my_stat_1d('./data/dwdy2m1d.dat',nb_point);
[y,dwdz2] = read_my_stat_1d('./data/dwdz2m1d.dat',nb_point);

diss = dudx2 + dudy2 + dudz2 + dvdx2 + dvdy2 + dvdz2 + dwdx2 + dwdy2 + dwdz2;
dissm = zeros(ym);
dissp = zeros(yp);
for i=1:ny
  dissm(i) = (1+mu_tm(i)/mu) * diss(i) / (re*utau_m*utau_m * re*utau_m*utau_m);
  dissp(i) = (1+mu_tp(i)/mu) * diss(nb_point-i+1) / (re*utau_p*utau_p * re*utau_p*utau_p);
end

[y,mudiss] = read_my_stat_1d('./data/mutudissm1d.dat',nb_point);
mudiss = mudiss - mu_t .* dudx2 - mu_t .* dvdx2 - mu_t .* dwdx2;
mudiss = mudiss - mu_t .* dudy2 - mu_t .* dvdy2 - mu_t .* dwdy2;
mudiss = mudiss - mu_t .* dudz2 - mu_t .* dvdz2 - mu_t .* dwdz2;
mudiss = mudiss - 2*dudx.*mududx - 2*dvdx.*mudvdx - 2*dwdx.*mudwdx;
mudiss = mudiss - 2*dudy.*mududy - 2*dvdy.*mudvdy - 2*dwdy.*mudwdy;
mudiss = mudiss - 2*dudz.*mududz - 2*dvdz.*mudvdz - 2*dwdz.*mudwdz;
mudiss = mudiss + mu_t.*dudx.*dudx + mu_t.*dvdx.*dvdx + mu_t.*dwdx.*dwdx;
mudiss = mudiss + mu_t.*dudy.*dudy + mu_t.*dvdy.*dvdy + mu_t.*dwdy.*dwdy;
mudiss = mudiss + mu_t.*dudz.*dudz + mu_t.*dvdz.*dvdz + mu_t.*dwdz.*dwdz;
mudiss = mudiss / mu;
mudissm = zeros(ym);
mudissp = zeros(yp);
for i=1:ny
  mudissm(i) = mudiss(i) / (re*utau_m*utau_m * re*utau_m*utau_m);
  mudissp(i) = mudiss(nb_point-i+1) / (re*utau_p*utau_p * re*utau_p*utau_p);
end

[y,ume]      = read_my_stat_1d('./data/um1d.dat',nb_point);
write_csv([ym,dissm,mudissm], "./csv/diss1.csv"," ");
write_csv([yp,dissp,mudissp], "./csv/diss2.csv"," ");

dy = zeros(ny,1);
dy(1)=2*y(1);
for i=2:ny
  dy(i)=2*(y(i)-sum(dy));
end
vcel2 = ((dx*dy(1:ny)*dz).^(2/3));
smod = sqrt(2).*sqrt(dudx2+dvdy2+dwdz2+0.5*(dudy.*dudy+dudy2+dudz2+dvdx2+dvdz2+dwdx2+dwdy2));
cs = mu_t(1:ny)./( 4*vcel2(1:ny).*smod(1:ny) );
write_csv([ym,um,mu_tm/mu,sqrt(cs)], "./csv/moy1.csv"," ");

exit;
