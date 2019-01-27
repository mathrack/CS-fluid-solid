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

#Graph settings
axis([5,150,0.,0.12])
xscale('log')
yscale('linear')
xlabel(r"$y^+$")
ylabel(r"$\varepsilon_\theta$")

# Fluid-solid interface
plot(array([0,0]),array([0.00000001,1000]),'--',color='grey')

# DNS, isoQ
refn = xlrd.open_workbook('../data/neuma.xls')
sh2 = refn.sheet_by_name('budget_tt')
y_ref = purge(sh2.col_values(0)[3:101])
tt_ref=-purge(sh2.col_values(1)[3:101])
plot(y_ref,tt_ref,'-',color='k',label=r'$isoQ$',markerfacecolor='none',markeredgecolor='k')
# LES, isoQ
lesn = xlrd.open_workbook('../../LES_database/wale_150_pr071/wale_150_pr071.xlsm')
les = lesn.sheet_by_name('scalar')
y_sol_les = purge(les.col_values(0)[1:33])
y_les = purge(les.col_values(0)[33:64])
tt_les = purge(les.col_values(6)[97:128])
plot(y_les, tt_les,'-o',color='k',markerfacecolor='none',markeredgecolor='k')

# DNS, G=1/2, G2=1, K=1/(G2*sqrt(G))=sqrt(2)
cht =  xlrd.open_workbook('../data/g05a1.xls')
sh3 = cht.sheet_by_name('budget_tt')
y_cht=purge(sh3.col_values(0)[3:101])
y_sol=purge(sh3.col_values(8)[3:135])
tt_cht=-purge(sh3.col_values(1)[3:101])
somm_conj=array(purge(sh3.col_values(5)[3:101]))
somm_conj[0:6]=0.
tt_cht=tt_cht+somm_conj
tt_sol=purge(sh3.col_values(10)[3:135])
plot(y_cht,tt_cht,'--',color='g',label=r'$G=\frac{1}{2}$')
#plot(y_sol,-tt_sol,'--',color='g')
tt_les = purge(les.col_values(6)[193:225])
#plot(y_sol_les,tt_les,'--+',color='g',markerfacecolor='none',markeredgecolor='g')
tt_les = purge(les.col_values(6)[225:256])
plot(y_les,tt_les,'--+',color='g',markerfacecolor='none',markeredgecolor='g')

# DNS, G=2, G2=1/2, K=1/(G2*sqrt(G))=sqrt(2)
cht =  xlrd.open_workbook('../data/g2a2.xls')
sh3 = cht.sheet_by_name('budget_tt')
y_cht=purge(sh3.col_values(0)[3:101])
y_sol=purge(sh3.col_values(8)[3:135])
tt_cht=-purge(sh3.col_values(1)[3:101])
somm_conj=array(purge(sh3.col_values(5)[3:101]))
somm_conj[0:6]=0.
tt_cht=tt_cht+somm_conj
tt_sol=purge(sh3.col_values(10)[3:135])
plot(y_cht,tt_cht,'-.',color='r',label=r'$G=2$')
#plot(y_sol,-tt_sol,'-.',color='r')
tt_les = purge(les.col_values(6)[641:673])
#plot(y_sol_les,tt_les,'-.x',color='r',markerfacecolor='none',markeredgecolor='r')
tt_les = purge(les.col_values(6)[673:704])
plot(y_les,tt_les,'-.x',color='r',markerfacecolor='none',markeredgecolor='r')

# DNS, isoT
refd = xlrd.open_workbook('../data/diric.xls')
sh2 = refd.sheet_by_name('budget_tt')
y_refd=purge(sh2.col_values(0)[3:101])
tt_refd=-purge(sh2.col_values(1)[3:101])
plot(y_refd,tt_refd,':',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')
# LES, isoT
tt_les = purge(les.col_values(6)[33:64])
plot(y_les, tt_les,':o',color='k',markerfacecolor='none',markeredgecolor='k')

legend(bbox_to_anchor=(1.00,1.0),numpoints=1)

savefig("all_diss_K2_fluid.png",bbox_inches='tight')
savefig("all_diss_K2_fluid.pdf",bbox_inches='tight')


