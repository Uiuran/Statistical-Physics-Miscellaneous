from importall import *
from matplotlib import rc

rc('text',usetex = True);
rc('font',family = 'serif');
alpha = 1.4 # parametro  'a' da lorentziana
x = np.linspace(-8.5,8.5,200000); # dominio da funcao
a = alpha*np.tan(np.pi*(np.random.uniform(0,1,1000000) - 0.5)); # transformacao de variavel aleatoria uniformemente distribuida em [0,1] para variavel exponencialmente distribuida

def histogram(data,bins, minimum = -7.5, maximum = 7.5):
   
    binsize = (maximum - minimum)/float(bins);
    binxpos = [minimum + binsize*j for j in np.arange(bins)];
    binypos = [];
    data.sort();
    for i in range(bins):
        s = data.searchsorted([minimum + i*binsize,minimum + (i+1)*binsize]);
        binguys = data[s[0]:s[1]];
        binypos.append(binguys.size);
    norm = float(binsize*data.size) #float(sum(binypos));
    binypos = [j/norm for j in binypos];
    plt.bar(binxpos,binypos,binsize,color='red');

plt.plot(x, alpha/(np.pi*(alpha**2 + x**2 )) );
histogram(a, 100)
plt.ylabel(r'$\rho(x)$', size = 20);
plt.legend((r' $a/(\pi(a^{2} + x^{2}))$',r'histograma 100 bins $10^{6}$ n\' umeros sorteados'));
plt.show()
