##################################################################
# Flageul C., Tiselj I., Benhamadouche S., Ferrand M.
# https://doi.org/10.1007/s10494-019-00008-0
##################################################################
from numpy import *
from pylab import *
import xlrd
matplotlib.rc('figure', figsize=(5.83,4.13))
matplotlib.rc('text', usetex = True)
size=16
size_legend=14
size_label=20
linewidth=1.5
markersize=10
matplotlib.rc('lines', linewidth=linewidth,markersize=markersize)
matplotlib.rc('font', size=size)
matplotlib.rc('axes', labelsize=size_label, titlesize=size)
matplotlib.rc('legend', fontsize=size_legend)

#
def purge(L):
    return array([x for x in L if x <> ''])

class myscalar(object):
     def __init__(self, wksht, ligne, ny):
         self.ny = ny
         self.name = wksht.col_values(0)[ligne]
         self.y      = purge(wksht.col_values(0)[1+ligne:1+ligne+self.ny-1])
         self.t      = purge(wksht.col_values(1)[1+ligne:1+ligne+self.ny-1])
         self.tt     = purge(wksht.col_values(2)[1+ligne:1+ligne+self.ny-1])
         self.ut     = purge(wksht.col_values(3)[1+ligne:1+ligne+self.ny-1])
         self.vt     = purge(wksht.col_values(4)[1+ligne:1+ligne+self.ny-1])
         self.wt     = purge(wksht.col_values(5)[1+ligne:1+ligne+self.ny-1])
         self.rdiss  = purge(wksht.col_values(6)[1+ligne:1+ligne+self.ny-1])
         self.ydiss  = purge(wksht.col_values(7)[1+ligne:1+ligne+self.ny-1])
         self.mudiss = purge(wksht.col_values(8)[1+ligne:1+ligne+self.ny-1])
         self.diss = - self.rdiss - self.mudiss

class case(object):
     def __init__(self, xlsfile, ncell, nsol, nscalar, pr):
         self.file = xlsfile
         self.ny = ncell+1
         self.nsol = nsol
         self.nscal = nscalar
         self.pr = pr
         self.wkbk = xlrd.open_workbook(self.file)
         self.Usht = self.wkbk.sheet_by_name('<U>')
         self.y = purge(self.Usht.col_values(0)[1:self.ny])
         self.u = purge(self.Usht.col_values(1)[1:self.ny])
         self.mu = purge(self.Usht.col_values(2)[1:self.ny])
         #self.cs = purge(self.Usht.col_values(3)[1:self.ny])
         self.Rijsht = self.wkbk.sheet_by_name('Rij')
         self.rxx = purge(self.Rijsht.col_values(1)[1:self.ny])
         self.ryy = purge(self.Rijsht.col_values(2)[1:self.ny])
         self.rzz = purge(self.Rijsht.col_values(3)[1:self.ny])
         self.rxy = purge(self.Rijsht.col_values(4)[1:self.ny])
         self.Disssht = self.wkbk.sheet_by_name('Diss.')
         self.rdiss = purge(self.Disssht.col_values(1)[1:self.ny])
         self.mudiss = purge(self.Disssht.col_values(2)[1:self.ny])
         self.diss = - self.rdiss - self.mudiss
         if self.nscal > 0:
             self.scalsht = self.wkbk.sheet_by_name('scalar')
             ligne=0
             self.scal=[myscalar(self.scalsht, i*(self.nsol+self.ny), (self.nsol+self.ny)) for i in range(self.nscal)]
         if self.nscal > 48:
             self.ttfs=zeros(49)
             for i in range(7):
                 for j in range(7):
                     self.ttfs[i+7*j] = self.scalsht.cell_value(42+j,16+i)
             self.epf_fs=zeros(49)
             for i in range(7):
                 for j in range(7):
                     self.epf_fs[i+7*j] = self.scalsht.cell_value(63+j,25+i)
             self.eps_fs=zeros(49)
             for i in range(7):
                 for j in range(7):
                     self.eps_fs[i+7*j] = self.scalsht.cell_value(70+j,25+i)
         else:
             self.ttfs=zeros(3)
             for i in range(3):
                 self.ttfs[i] = self.scalsht.cell_value(30+i,16)
             self.epf_fs=zeros(3)
             for i in range(3):
                 self.epf_fs[i] = (self.scalsht.cell_value(39+i,16)+self.scalsht.cell_value(42+i,16))/self.pr
             self.eps_fs=zeros(3)
             for i in range(3):
                 self.eps_fs[i] = (self.scalsht.cell_value(39+i,16)+self.scalsht.cell_value(45+i,16))/(self.pr*self.scalsht.cell_value(21+i,16))

# LES results
wale_395_071 = case('../../LES_database/wale_395_pr071/wale_395_pr071.xlsm', 42, 31, 51, 0.71)
wale_395_1 = case('../../LES_database/wale_395_pr1/wale_395_pr1.xlsm', 42, 31, 5, 1.)

wale_1020_071 = case('../../LES_database/wale_1020_pr071/wale_1020_pr071.xlsm', 86, 31, 5, 0.71)
wale_1020_1 = case('../../LES_database/wale_1020_pr1/wale_1020_pr1.xlsm', 86, 31, 5, 1.)

#Graph settings
xscale('linear')
yscale('log')
axis([-10,10,0.001,0.2])
xlabel(r"$y^+$")
ylabel(r"$\varepsilon_\theta$")
title(r"Cases $G=1$ and $K=1$")

# Fluid-solid interface
plot(array([0,0]),array([0.00000001,1000]),'--',color='grey')

iscal = 2+3*7+3
plot(append(wale_395_071.scal[iscal].y[:31],0),append(-wale_395_071.scal[iscal].diss[:31],wale_395_071.eps_fs[iscal-2]),'-',color='k',label=r'$Re_\tau|Pr=395|0.71$')
plot(insert(wale_395_071.scal[iscal].y[31:],0,0),insert(-wale_395_071.scal[iscal].diss[31:],0,wale_395_071.epf_fs[iscal-2]),'-',color='k')

iscal = 2
plot(append(wale_395_1.scal[iscal].y[:31],0),append(-wale_395_1.scal[iscal].diss[:31],wale_395_1.eps_fs[iscal-2]),'--',color='r',label=r'$Re_\tau|Pr=395|1$')
plot(insert(wale_395_1.scal[iscal].y[31:],0,0),insert(-wale_395_1.scal[iscal].diss[31:],0,wale_395_1.epf_fs[iscal-2]),'--',color='r')

iscal = 2
plot(append(wale_1020_071.scal[iscal].y[:31],0),append(-wale_1020_071.scal[iscal].diss[:31],wale_1020_071.eps_fs[iscal-2]),'-.',color='g',label=r'$Re_\tau|Pr=1020|0.71$')
plot(insert(wale_1020_071.scal[iscal].y[31:],0,0),insert(-wale_1020_071.scal[iscal].diss[31:],0,wale_1020_071.epf_fs[iscal-2]),'-.',color='g')

iscal = 2
plot(append(wale_1020_1.scal[iscal].y[:31],0),append(-wale_1020_1.scal[iscal].diss[:31],wale_1020_1.eps_fs[iscal-2]),':',color='b',label=r'$Re_\tau|Pr=1020|1$')
plot(insert(wale_1020_1.scal[iscal].y[31:],0,0),insert(-wale_1020_1.scal[iscal].diss[31:],0,wale_1020_1.epf_fs[iscal-2]),':',color='b')

legend(bbox_to_anchor=(0.45,0.5),numpoints=1)

savefig("diss_g1k1.png",bbox_inches='tight')
savefig("diss_g1k1.pdf",bbox_inches='tight')
