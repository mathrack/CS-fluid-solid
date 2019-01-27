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

class case(object):

     def __init__(self, xlsfile, ncell):
         self.file = xlsfile
         self.ny = ncell
         self.wkbk = xlrd.open_workbook(self.file)
         self.Usht = self.wkbk.sheet_by_name('<U>')
         self.y = purge(self.Usht.col_values(0)[1:self.ny])
         self.u = purge(self.Usht.col_values(1)[1:self.ny])
         self.mu = purge(self.Usht.col_values(2)[1:self.ny])
         self.Rijsht = self.wkbk.sheet_by_name('Rij')
         self.rxx = purge(self.Rijsht.col_values(1)[1:self.ny])
         self.ryy = purge(self.Rijsht.col_values(2)[1:self.ny])
         self.rzz = purge(self.Rijsht.col_values(3)[1:self.ny])
         self.rxy = purge(self.Rijsht.col_values(4)[1:self.ny])
         self.Disssht = self.wkbk.sheet_by_name('Diss.')
         self.rdiss = purge(self.Disssht.col_values(1)[1:self.ny])
         self.mudiss = purge(self.Disssht.col_values(2)[1:self.ny])

# DNS results
tmp_nomod = xlrd.open_workbook('../1020_blaz_wale_foam/1020_blaz_wale_foam.xlsx')
# DNS
tmp_dns = tmp_nomod.sheet_by_name('DNS')
y_dns = purge(tmp_dns.col_values(1)[1:224])
u_dns = purge(tmp_dns.col_values(2)[1:224])
rxx_dns = purge(tmp_dns.col_values(3)[1:224])
rzz_dns = purge(tmp_dns.col_values(4)[1:224])
ryy_dns = purge(tmp_dns.col_values(2)[227:450])
rxy_dns = purge(tmp_dns.col_values(3)[227:450])
diss_dns = purge(tmp_dns.col_values(3)[458:681])
# Close the excel file
tmp_nomod.release_resources()

nomod = case('../1020_ygrid_nomod/1020_ygrid_nomod.xlsx', 86)
smag = case('../1020_ygrid_smag/1020_ygrid_smag.xlsx', 86)
wale = case('../1020_ygrid_wale/1020_ygrid_wale.xlsx', 86)
dyn0 = case('../1020_ygrid_dyn0/1020_ygrid_dyn0.xlsx', 86)
dyn1 = case('../1020_ygrid_dyn1/1020_ygrid_dyn1.xlsx', 86)
dyn0b = case('../1020_ygrid_dyn0b2/1020_ygrid_dyn0b2.xlsx', 86)
dyn1b = case('../1020_ygrid_dyn1b2/1020_ygrid_dyn1b2.xlsx', 86)

#Graph settings
xscale('log')
yscale('linear')
axis([0.4,1051,0.,1.])
xlabel(r"$y^+$")
ylabel(r"$R_{xy}$")

plot(y_dns,rxy_dns,'-',color='k',label=r'$DNS$')
plot(nomod.y,nomod.rxy,'--',color='r',label=r'$No.mod.$')
plot(smag.y,smag.rxy,'--',color='g',label=r'$Smag.$')
plot(wale.y,wale.rxy,'--',color='b',label=r'$Wale$')
plot(dyn0.y,dyn0.rxy,'-.',color='r',label=r'$Dyn.(def.)$')
plot(dyn1.y,dyn1.rxy,'-.',color='b',label=r'$Dyn.(clip)$')
plot(dyn0b.y,dyn0b.rxy,':',color='r',label=r'$Dyn.(def.)+<>$')
plot(dyn1b.y,dyn1b.rxy,':',color='b',label=r'$Dyn.(clip)+<>$')

legend(bbox_to_anchor=(1.6,0.9),numpoints=1)

savefig("rxy.png",bbox_inches='tight')
savefig("rxy.pdf",bbox_inches='tight')
