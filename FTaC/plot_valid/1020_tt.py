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

def purge(L):
	return array([x for x in L if x <> ''])

#Graph settings
xscale('linear')
yscale('linear')
axis([0.,10,0.,8.])
xlabel(r"$y^+$")
ylabel("$\overline{T'^2}$")

# Fluid-solid interface
plot(array([0,0]),array([0.00000001,1000]),'--',color='grey')

# DNS, isoT
dnsn = xlrd.open_workbook('../data/1020.xlsx')
dns = dnsn.sheet_by_name('scalar')
y_ref=purge(dns.col_values(1)[165:389])
tt_ref=purge(dns.col_values(3)[165:389])
plot(y_ref,tt_ref*tt_ref,':',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')

# LES, isoT
lesn = xlrd.open_workbook('../../LES_database/wale_1020_pr071/wale_1020_pr071.xlsm')
les = lesn.sheet_by_name('scalar')
y_les = purge(les.col_values(0)[32:118])
tt_les = purge(les.col_values(2)[32:118])
plot(y_les, tt_les,':o',color='k',markerfacecolor='none',markeredgecolor='k')

legend(bbox_to_anchor=(0.45,0.7),numpoints=1)
savefig("1020_tt.png",bbox_inches='tight')
savefig("1020_tt.pdf",bbox_inches='tight')

xscale('log')
axis([5.,1020.,0.,8.])
savefig("1020_tt_xlog.png",bbox_inches='tight')
savefig("1020_tt_xlog.pdf",bbox_inches='tight')
