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
	return [x for x in L if x <> '']

#Graph settings
xscale('linear')
yscale('linear')
axis([-10,10,0.,9])
xlabel(r"$y^+$")
ylabel("$\overline{T'^2}$")

# Fluid-solid interface
plot(array([0,0]),array([0.00000001,1000]),'--',color='grey')

# DNS, isoQ
refn = xlrd.open_workbook('../data/neuma.xls')
sh2 = refn.sheet_by_name('fluctuations')
y_ref=purge(sh2.col_values(9)[5:101])
tt_ref=purge(sh2.col_values(16)[5:101])
plot(y_ref,tt_ref,'-',color='k',label=r'$isoQ$',markerfacecolor='none',markeredgecolor='k')
# LES
lesn = xlrd.open_workbook('../../LES_database/wale_150_pr071/wale_150_pr071.xlsm')
les = lesn.sheet_by_name('scalar')
y_sol_les = purge(les.col_values(0)[1:64])
y_les = purge(les.col_values(0)[33:64])
tt_les = purge(les.col_values(2)[97:128])
plot(y_les, tt_les,'-o',color='k',markerfacecolor='none',markeredgecolor='k')

# DNS, G=1/2, G2=1, K=1/(G2*sqrt(G))=sqrt(2)
cht =  xlrd.open_workbook('../data/g05a1.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'--',color='g',label=r'$G=\frac{1}{2}$')
plot(y_sol,tt_sol,'--',color='g')
# LES
tt_les = purge(les.col_values(2)[193:256])
plot(y_sol_les,tt_les,'--+',color='g',markerfacecolor='none',markeredgecolor='g')

# DNS, G=2, G2=1/2, K=1/(G2*sqrt(G))=sqrt(2)
cht =  xlrd.open_workbook('../data/g2a2.xls')
sh3 = cht.sheet_by_name('fluctuations')
y_cht=purge(sh3.col_values(9)[5:101])
y_sol=purge(sh3.col_values(28)[5:135])
tt_cht=purge(sh3.col_values(16)[5:101])
tt_sol=purge(sh3.col_values(29)[5:135])
plot(y_cht,tt_cht,'-.',color='r',label=r'$G=2$')
plot(y_sol,tt_sol,'-.',color='r')
# LES
tt_les = purge(les.col_values(2)[641:704])
plot(y_sol_les,tt_les,'-.x',color='r',markerfacecolor='none',markeredgecolor='r')

# DNS, isoT
refd = xlrd.open_workbook('../data/diric.xls')
sh2d= refd.sheet_by_name('fluctuations')
y_refd=purge(sh2d.col_values(9)[5:101])
tt_refd=purge(sh2d.col_values(16)[5:101])
plot(y_refd,tt_refd,':',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')
# LES, isoT
tt_les = purge(les.col_values(2)[33:64])
plot(y_les, tt_les,':o',color='k',markerfacecolor='none',markeredgecolor='k')

text(-3, 7.5, r'$Solid$', horizontalalignment='center',verticalalignment='center',color='k')
text(2, 7.5, r'$Fluid$', horizontalalignment='center',verticalalignment='center',color='k')

legend(bbox_to_anchor=(0.4,0.7),numpoints=1)
savefig("all_tt_K2.png",bbox_inches='tight')
savefig("all_tt_K2.pdf",bbox_inches='tight')
