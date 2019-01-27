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
xscale('log')
yscale('linear')
axis([5,400,0.,9.])
xlabel(r"$y^+$")
ylabel("$\overline{T'^2}$")

# Fluid-solid interface
plot(array([0,0]),array([0.00000001,1000]),'--',color='grey')

# DNS, isoQ
dnsn = xlrd.open_workbook('../data/395_dns.xlsx')
dns = dnsn.sheet_by_name('fluctuations')
y_ref=purge(dns.col_values(29)[2:131])
tt_ref=purge(dns.col_values(53)[2:131])
plot(y_ref,tt_ref,'-',color='k',label=r'$isoQ$',markerfacecolor='none',markeredgecolor='k')
# LES
lesn = xlrd.open_workbook('../data/395_les.xlsx')
les = lesn.sheet_by_name('scalar')
y_sol_les = purge(les.col_values(0)[1:85])
y_les = purge(les.col_values(0)[43:85])
tt_les = purge(les.col_values(2)[128:170])
plot(y_les, tt_les,'-o',color='k',markerfacecolor='none',markeredgecolor='k')

# DNS, G=1, G2=1, K=1/(G2*sqrt(G))=1
y_cht=y_ref
tt_cht=purge(dns.col_values(56)[2:131])
y_sol=purge(dns.col_values(59)[2:259])
tt_sol=purge(dns.col_values(60)[2:259])
plot(y_cht,tt_cht,'--',color='g',label=r'$G=G_2=1$')
plot(y_sol,tt_sol,'--',color='g')
# LES
tt_les = purge(les.col_values(2)[171:255])
plot(y_sol_les,tt_les,'--+',color='g',markerfacecolor='none',markeredgecolor='g')

# DNS, isoT
y_refd=y_ref
tt_refd=purge(dns.col_values(52)[2:131])
plot(y_refd,tt_refd,':',color='k',label=r'$isoT$',markerfacecolor='none',markeredgecolor='k')
# LES, isoT
tt_les = purge(les.col_values(2)[43:85])
plot(y_les, tt_les,':o',color='k',markerfacecolor='none',markeredgecolor='k')

legend(bbox_to_anchor=(0.55,0.6),numpoints=1)
savefig("395_tt_fluid.png",bbox_inches='tight')
savefig("395_tt_fluid.pdf",bbox_inches='tight')
