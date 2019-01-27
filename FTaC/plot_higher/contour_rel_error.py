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
     def __init__(self, xlsfile, ncell, nsol, nscalar):
         self.file = xlsfile
         self.ny = ncell+1
         self.nsol = nsol
         self.nscal = nscalar
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
         self.gg=zeros([7,7])
         for i in range(7):
             for j in range(7):
                 self.gg[i][j] = self.scalsht.cell_value(21+j,16+i)
         self.kk=zeros([7,7])
         for i in range(7):
             for j in range(7):
                 self.kk[i][j] = self.scalsht.cell_value(28+j,16+i)
         self.ttfs=zeros(49)
         for i in range(7):
             for j in range(7):
                 self.ttfs[i+7*j] = self.scalsht.cell_value(42+j,16+i)
         self.epf_fs=zeros(49)
         self.epf_fs_mat=zeros([7,7])
         for i in range(7):
             for j in range(7):
                 self.epf_fs[i+7*j] = self.scalsht.cell_value(63+j,25+i)
                 self.epf_fs_mat[i][j] = self.epf_fs[i+7*j]
         self.eps_fs=zeros(49)
         self.eps_fs_mat=zeros([7,7])
         for i in range(7):
             for j in range(7):
                 self.eps_fs[i+7*j] = self.scalsht.cell_value(70+j,25+i)
                 self.eps_fs_mat[i][j] = self.eps_fs[i+7*j]

# LES results
wale = case('../../LES_database/wale_395_pr071/wale_395_pr071.xlsm', 42, 31, 51)

#Graph settings
xscale('linear')
yscale('linear')
xlabel(r"$log_{10}(G)$")
ylabel(r"$log_{10}(K)$")
title(r"Contour plot of the relative error, in \%")

fit = (0.00197363782954315/0.0247008006107844) * (wale.gg**0.225405802) * (wale.kk**1.901499052)
fit = 1/wale.gg + divide(wale.kk**2 - 1/wale.gg, 1 +  fit)

rel_error = 1.0 - divide(fit, divide(wale.eps_fs_mat, wale.epf_fs_mat) )

levels = [-15,-10,-4,-1,1,4,10,15]
CS = contour(log10(wale.gg),log10(wale.kk),100*rel_error,levels,colors='k')
fmt = r'%d'
clabel(CS, inline=True, fmt=fmt, fontsize=20)

savefig("contour_rel_error.png",bbox_inches='tight')
savefig("contour_rel_error.pdf",bbox_inches='tight')
