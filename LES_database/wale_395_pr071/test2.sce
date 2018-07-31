nb_point=146;
nb_sol=31;
nx=81;
nz=81;
lx=6.283;
ly=2.;
lz=3.142;
mu=2.0/14124.0;
re=1.0/mu;
pr=0.71;
dx = lx / nx;
dz = lz / nz;

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

// This will be used to compute cos(alpha)
yc = y;
ypp = zeros(nb_point+1,1);
dy = zeros(nb_point,1);
ypp(1) = -150./400.;
for i = 1:nb_point
  dy(i) = 2.*(yc(i)-ypp(i));
  ypp(i+1) = ypp(i) + dy(i);
end
alpha = zeros(nb_point-1,1);
for i = 1:nb_point-1
  alpha(i) = (yc(i+1)-ypp(i+1)) / (yc(i+1)-yc(i));
end
yi = yc(2:nb_point-1);
function [dery,var_p] = compute_grad(var,nn,alpha,dy,g2)
  dery = zeros(nn-2,1);
  var_p = zeros(nn-1);
  for jj = 1:30
    var_p(jj) = ( var(jj)/(1.-alpha(jj)) + var(jj+1)/alpha(jj) ) / ( 1./(1.-alpha(jj)) + 1./alpha(jj) );
  end
  jj = 31; var_p(jj) = ( var(jj)/(1.-alpha(jj)) + (1./g2)*var(jj+1)/alpha(jj) ) / ( 1./(1.-alpha(jj)) + (1./g2)/alpha(jj) );
  for jj = 32:114
    var_p(jj) = ( var(jj)/(1.-alpha(jj)) + var(jj+1)/alpha(jj) ) / ( 1./(1.-alpha(jj)) + 1./alpha(jj) );
  end
  jj = 115; var_p(jj) = ( (1./g2)*var(jj)/(1.-alpha(jj)) + var(jj+1)/alpha(jj) ) / ( (1./g2)/(1.-alpha(jj)) + 1./alpha(jj) );
  for jj = 116:nn-1
    var_p(jj) = ( var(jj)/(1.-alpha(jj)) + var(jj+1)/alpha(jj) ) / ( 1./(1.-alpha(jj)) + 1./alpha(jj) );
  end
  dery(1:nn-2) = (var_p(2:nn-1) - var_p(1:nn-2) ) ./ dy(2:nn-1);
endfunction;
function outp = compute_integral(var,dy)
  outp = 0.;
  for jj = 1:31
    outp = outp + var(jj)*dy(jj);
  end
endfunction;

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

utau_m=sqrt( (abs(dudy(nb_sol+1)))/re );
utau_p=sqrt( (abs(dudy(nb_point-nb_sol)))/re );
utau=sqrt(0.5*(utau_m*utau_m+utau_p*utau_p));

// Print the estimation of the friction Reynolds number
re*utau
// Use the theoretical value (imposed pressure gradient)
utau = 395. / re;
utau_m = utau;
utau_p = utau;

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
  dissm(i) = (1 + mu_tm(i)/mu) .* diss(i) / (re*utau_m*utau_m * re*utau_m*utau_m);
  dissp(i) = (1 + mu_tp(i)/mu) .* diss(nb_point-i+1) / (re*utau_p*utau_p * re*utau_p*utau_p);
end

[y,mudiss] = read_my_stat_1d('./data/mutudissm1d.dat',nb_point);
mudiss = mudiss - mu_t.*diss;
mudiss = mudiss - 2*(dudx.*mududx + dvdx.*mudvdx + dwdx.*mudwdx);
mudiss = mudiss - 2*(dudy.*mududy + dvdy.*mudvdy + dwdy.*mudwdy);
mudiss = mudiss - 2*(dudz.*mududz + dvdz.*mudvdz + dwdz.*mudwdz);
mudiss = mudiss + mu_t.*(dudx.*dudx + dudy.*dudy + dudz.*dudz + dvdx.*dvdx + dvdy.*dvdy + dvdz.*dvdz + dwdx.*dwdx + dwdy.*dwdy + dwdz.*dwdz);
mudiss = mudiss / mu;
mudissm = zeros(ym);
mudissp = zeros(yp);
for i=1:ny
  mudissm(i) = mudiss(i) / (re*utau_m*utau_m * re*utau_m*utau_m);
  mudissp(i) = mudiss(nb_point-i+1) / (re*utau_p*utau_p * re*utau_p*utau_p);
end

// Diric + Neuma + Robin + Conjg
nb_head = 51;
head=["diric","neuma","11","12","13","14","15","16","17","21","22","23","24","25","26","27","31","32","33","34","35","36","37","41","42","43","44","45","46","47","51","52","53","54","55","56","57","61","62","63","64","65","66","67","71","72","73","74","75","76","77"];
gggg=[ 1.    , 1.    , 0.1,0.1,0.1,0.1,0.1,0.1,0.1, 0.2,0.2,0.2,0.2,0.2,0.2,0.2, 0.5,0.5,0.5,0.5,0.5,0.5,0.5, 1.,1.,1.,1.,1.,1.,1., 2.,2.,2.,2.,2.,2.,2., 5.,5.,5.,5.,5.,5.,5., 10.,10.,10.,10.,10.,10.,10.];
kkkk=[ 1.    , 1.    , 0.1,0.2,0.5,1.,2.,5.,10., 0.1,0.2,0.5,1.,2.,5.,10., 0.1,0.2,0.5,1.,2.,5.,10., 0.1,0.2,0.5,1.,2.,5.,10., 0.1,0.2,0.5,1.,2.,5.,10., 0.1,0.2,0.5,1.,2.,5.,10., 0.1,0.2,0.5,1.,2.,5.,10.];

myii = zeros(51,1);
myjj = zeros(51,1);
mygg = zeros(7,7);
mykk = zeros(7,7);
myg2 = zeros(7,7);
mytt = zeros(7,7);
myepsf = zeros(7,7);
myepss = zeros(7,7);
mydyttf = zeros(7,7);
mydytts = zeros(7,7);
mydtdtxz = zeros(7,7);
mydtdtyf = zeros(7,7);
mydtdtys = zeros(7,7);
myerror = zeros(7,7);
myepss_sur_epsf = zeros(7,7);

// mean + variance + flux + grad + diss + mut*grad + mut*diss
nb_tail = 1+1+3+3+6+3+1;
tail=["_mean.dat",..
      "_var.dat",..
      "_ux.dat","_uy.dat","_uz.dat",..
      "_dtdx.dat","_dtdy.dat","_dtdz.dat",..
      "_dtdx2.dat","_dtdy2.dat","_dtdz2.dat","_dtdxy.dat","_dtdxz.dat","_dtdyz.dat",..
      "_mu_dtdx.dat","_mu_dtdy.dat","_mu_dtdz.dat",..
      "_mu_diss.dat"];
// Read and write
for i = 1:nb_head
  j=1;
  [y,t_mean] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_var]  = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_ux]   = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_uy]   = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_uz]   = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdx] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdy] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdz] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdx2] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdy2] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdz2] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdxy] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdxz] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_dtdyz] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_mudtdx] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_mudtdy] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_mudtdz] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;
  [y,t_mudiss] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point); j=j+1;

  t_ux = t_ux - t_mean.*ume;
  t_uy = t_uy - t_mean.*vme;
  t_uz = t_uz - t_mean.*wme;
  t_mudiss = t_mudiss - mu_t .* (t_dtdx2+t_dtdy2+t_dtdz2) / 0.5;
  t_mudiss = t_mudiss - 2*(t_dtdx.*t_mudtdx + t_dtdy.*t_mudtdy + t_dtdz.*t_mudtdz);
  t_mudiss = t_mudiss + mu_t .* (t_dtdx.*t_dtdx + t_dtdy.*t_dtdy + t_dtdz.*t_dtdz) / 0.5;

     t_meanm = zeros(ny,1);
      t_varm = zeros(ny,1);
       t_uxm = zeros(ny,1);
       t_uym = zeros(ny,1);
       t_uzm = zeros(ny,1);
     t_dissm = zeros(ny,1);
    t_ydissm = zeros(ny,1);
   t_mudissm = zeros(ny,1);
     t_meanp = zeros(ny,1);
      t_varp = zeros(ny,1);
       t_uxp = zeros(ny,1);
       t_uyp = zeros(ny,1);
       t_uzp = zeros(ny,1);
     t_dissp = zeros(ny,1);
    t_ydissp = zeros(ny,1);
   t_mudissp = zeros(ny,1);
  for k=1:ny
     t_meanm(k) = t_mean(k);
      t_varm(k) = t_var(k);
       t_uxm(k) = t_ux(k);
       t_uym(k) = t_uy(k);
       t_uzm(k) = t_uz(k);
     t_dissm(k) = (1 + mu_tm(k)/(0.5*mu))*(t_dtdx2(k)+t_dtdy2(k)+t_dtdz2(k));
    t_ydissm(k) = (1 + mu_tm(k)/(0.5*mu))*t_dtdy2(k);
   t_mudissm(k) = t_mudiss(k);
     t_meanp(k) = t_mean(nb_point-k+1);
      t_varp(k) = t_var(nb_point-k+1);
       t_uxp(k) = t_ux(nb_point-k+1);
       t_uyp(k) = t_uy(nb_point-k+1);
       t_uzp(k) = t_uz(nb_point-k+1);
     t_dissp(k) = (1 + mu_tp(k)/(0.5*mu))*(t_dtdx2(nb_point-k+1)+t_dtdy2(nb_point-k+1)+t_dtdz2(nb_point-k+1));
    t_ydissp(k) = (1 + mu_tp(k)/(0.5*mu))*t_dtdy2(nb_point-k+1);
   t_mudissp(k) = t_mudiss(nb_point-k+1);
  end
  // This is the estimated friction temperature
  t_tau_m=(abs(t_dtdy(nb_sol+1)))/(re*pr*utau_m);
  t_tau_p=(abs(t_dtdy(nb_point-nb_sol)))/(re*pr*utau_p);
  // This is the theoretical friction temperature
  t_tau_m=1.0/(re*pr*utau_m);
  t_tau_p=1.0/(re*pr*utau_p);

// Rescale solid dissipation with thermal diffusivity ratio
t_dissm(1:nb_sol)=t_dissm(1:nb_sol)/gggg(i);
t_dissp(1:nb_sol)=t_dissp(1:nb_sol)/gggg(i);
t_ydissm(1:nb_sol)=t_ydissm(1:nb_sol)/gggg(i);                                                   
t_ydissp(1:nb_sol)=t_ydissp(1:nb_sol)/gggg(i);
t_mudissm(1:nb_sol)=t_mudissm(1:nb_sol)/gggg(i);
t_mudissp(1:nb_sol)=t_mudissp(1:nb_sol)/gggg(i);

  for k=1:ny
    if (t_dissm(k) == 0) then
      t_ydissm(k) = 0;
    else
      t_ydissm(k) = t_ydissm(k) / t_dissm(k);
    end
    if (t_dissp(k) == 0) then
      t_ydissp(k) = 0;
    else
      t_ydissp(k) = t_ydissp(k) / t_dissp(k);
    end
  end

  for k=1:ny
    if (t_dissm(k) == 0) then
      t_ydissm(k) = 0;
    else
      t_ydissm(k) = t_ydissm(k) * t_dissm(k);
    end
    if (t_dissp(k) == 0) then
      t_ydissp(k) = 0;
    else
      t_ydissp(k) = t_ydissp(k) * t_dissp(k);
    end
  end

  if ~(i==1) & ~(i==2) then
    nb_flu=nb_sol+1;
    g2 = 1/(kkkk(i)*sqrt(gggg(i)));
    //
      myii(i) = 1+floor((i-3)/7);
      myjj(i) = modulo(i-3,7)+1;
      mygg(myii(i),myjj(i)) = gggg(i);
      mykk(myii(i),myjj(i)) = kkkk(i);
      myg2(myii(i),myjj(i)) = g2;
    //
    t_tau=t_tau_m; u_tau=utau_m;
    ys=ym(nb_sol);
    yf=ym(nb_flu);
    // For temperature variance
    tts=t_varm(nb_sol)/(t_tau**2);
    ttf=t_varm(nb_flu)/(t_tau**2);
    eps=(t_dissm(nb_sol)+t_mudissm(nb_sol)/mu)/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    epf=(t_dissm(nb_flu)+t_mudissm(nb_flu)/mu)/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    mat=[1., -g2; yf, -ys];
    rhs1=2*pr*epf*yf-2*g2*gggg(i)*pr*eps*ys;
    rhs2=ttf+pr*epf*(yf**2)-tts-gggg(i)*pr*eps*(ys**2);
    aa=linsolve(mat,-[rhs1;rhs2]);
    // For wall-parallel dissipation
    xzeps=(t_dissm(nb_sol)-t_ydissm(nb_sol))/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    xzepf=(t_dissm(nb_flu)-t_ydissm(nb_flu))/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    mat=[1., -g2; yf, -ys];
    rhs1=0.;
    rhs2=pr*xzepf-gggg(i)*pr*xzeps;
    bb=linsolve(mat,-[rhs1;rhs2]);
    // Temperature variance at the fluid-solid interface
    af=aa(1);
    as=aa(2);
    //
    // Temperature variance : ttf-af*yf+pr*epf*(yf**2); OR tts-as*ys+pr*eps*(ys**2);
    //
      mytt(myii(i),myjj(i)) = ttf-af*yf+pr*epf*(yf**2);
    //

    //
    // Derivative of temperature variance
    //
    //   on the fluid side : af-2*pr*epf*yf; OR (as-2*gggg(i)*pr*eps*ys)*g2;
      mydyttf(myii(i),myjj(i)) = af-2*pr*epf*yf;
    //   on the solid side : (af-2*pr*epf*yf)/g2; OR as-2*gggg(i)*pr*eps*ys;
      mydytts(myii(i),myjj(i)) = as-2*gggg(i)*pr*eps*ys;

    //
    // Wall-parallel dissipation
    //
    af=bb(1);
    as=bb(2);
    //    on the fluid side : xzepf-af*yf/pr; OR gggg(i)*xzeps-as*ys/pr;
    //    on the solid side : xzepf/gggg(i)-af*yf/pr; OR xzeps-as*ys/pr;
    //
    // Reconstruction of dtdx*dtdx+dtdz*dtdz
    //
    xzeps=(t_dissm(nb_sol)-t_ydissm(nb_sol))/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    xzepf=(t_dissm(nb_flu)-t_ydissm(nb_flu))/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    //    on the fluid or solid side : pr*xzepf-bb(1)*yf; OR gggg(i)*pr*xzeps-bb(2)*ys;
    dtdxdtdx_dtdzdtdz = (pr*xzepf-bb(1)*yf);
      mydtdtxz(myii(i),myjj(i)) = dtdxdtdx_dtdzdtdz;

    //
    // Reconstruction of dtdy*dtdy
    //
    yeps=(t_ydissm(nb_sol))/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    yepf=(t_ydissm(nb_flu))/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    //   cos(alpha)
    cosalpha = 0.5*aa(2)/sqrt(tts*yeps*gggg(i)*pr);
    //    on the fluid side
    dtdydtdy_flu = ( (aa(1)-2*pr*epf*yf)/(2*cosalpha*sqrt(ttf-aa(1)*yf+pr*epf*(yf**2))) )**2;
      mydtdtyf(myii(i),myjj(i)) = dtdydtdy_flu;
    //    on the solid side
    dtdydtdy_sol = ( (aa(2)-2*gggg(i)*pr*eps*ys)/(2*cosalpha*sqrt(ttf-aa(1)*yf+pr*epf*(yf**2))) )**2;
      mydtdtys(myii(i),myjj(i)) = dtdydtdy_sol;

    //
    // Reconstruction of ratios and error
    //
    dtdt_flu = dtdydtdy_flu + dtdxdtdx_dtdzdtdz; myepsf(myii(i),myjj(i))=dtdt_flu/pr;
    dtdt_sol = dtdydtdy_sol + dtdxdtdx_dtdzdtdz; myepss(myii(i),myjj(i))=dtdt_sol/pr/gggg(i);
    // flu / sol
    flu_sur_sol = (kkkk(i)**2) * dtdydtdy_flu / dtdt_flu + (1 - dtdydtdy_flu / dtdt_flu) / gggg(i);
    // sol / flu
    sol_sur_flu = dtdydtdy_sol / dtdt_sol / (kkkk(i)**2) + (1 - dtdydtdy_sol / dtdt_sol) * gggg(i);
    // error
      myerror(myii(i),myjj(i)) = flu_sur_sol * sol_sur_flu;

    myepss_sur_epsf(myii(i),myjj(i)) = (dtdt_sol/pr/gggg(i)) / (dtdt_flu/pr);

//    [i, ttf-aa(1)*yf+pr*epf*(yf**2), dtdt_flu/pr, dtdt_sol/pr/gggg(i), myerror(myii(i),myjj(i))]

    s1 = mydytts(myii(i),myjj(i))/(2*gggg(i)*pr);
    s2 = 395*compute_integral(t_dissm,dy)/(pr*(re**2)*(t_tau**2)*(u_tau**2));
    [i, gggg(i), kkkk(i), s1, s2, 200*(s1-s2)/(s1+s2)]

  end

end

i=35
  j=2; [y,t_var]  = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point);
  j=10; [y,t_dtdy2] = read_my_stat_1d(strcat(["./data/",head(i),tail(j)]),nb_point);
  [dery,varp] = compute_grad(t_var, nb_point, alpha, dy, 1/(kkkk(i)*sqrt(gggg(i))));
  scf(); plot(yi,0.5*dery./(sqrt(t_var(2:nb_point-1)).*sqrt(t_dtdy2(2:nb_point-1))),'-o');
