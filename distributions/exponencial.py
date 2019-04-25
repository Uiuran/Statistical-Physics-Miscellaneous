from importall import *
from matplotlib import rc

rc('text',usetex = True);
rc('font',family = 'serif');
alpha = 1.4 # media da exponencial
x = np.linspace(0,10.0,100000); # dominio da funcao  exp{-x}
a = np.log(1.0/(np.random.uniform(0,1,1000)) )/alpha; # transformacao de variavel aleatoria uniformemente distribuida em [0,1] para variavel exponencialmente distribuida

plt.plot(x, alpha*np.exp(-alpha*x)); # plota forma exponencial
plt.hist(a,50, normed= True); # plota histograma normalizando pela quantidade de dados
plt.ylabel(r'$\rho(x)$', size = 20);
plt.legend((r' $e^{-x}$',r'histograma 50 bins $10^{3}$ n\' umeros sorteados'));
plt.show()
