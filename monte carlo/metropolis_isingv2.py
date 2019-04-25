import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
from matplotlib import rc

L = 100; # tamanho das dimensoes da grade
N = L*L; # numero de spins na grade 
J = 1.0 #coeficiente de spins cruzados
H = 0.0; # coeficiente de spin i
Time = 100000000; # tempo total
Meitempo = range(100000);
num = 10; # numero de realizacoes da simulacao de ising
#x =  np.linspace(0,9999,10000);

T = np.linspace(0.1,5.0,50) # temperaturas para simular dependencia da magnetizacao media termalizada (tempo total 200000 por realizacao, descarte das primeiras 100000)

S = np.random.rand(N); # condicoes iniciais
for i in range(len(S)):
    if S[i] > 0.5:
        S[i] = 1;
    else:
        S[i] = -1;
Sini = S.copy();

m = np.zeros((num,Time/2)); # magnetizacao
vm = np.zeros(Time/2); # variancia da magnetizacao, num realizacoes
mm = np.zeros(Time/2); # media da magnetizacao, num realizacoes

# Energia inicial
E = 0.0;
for j in range(len(S)):
    b = np.floor(j/L)*L;
    E = E -J*S[j]*(S[[np.mod(L+j,N),np.mod(j+1-b,L)]].sum());
E = E - H*S.sum();
Eini = E;
ef = np.zeros(num); # energias finais


## Simulacao ##
for temp in T:

    ## Loop para realizacoes em temperatura temp ##
    for i in range(num):        
        print temp,i
        for t in range(Time):

            a = np.random.randint(N);
            b = np.floor(a/L)*L;
            v = S[[np.mod(a-L,N) , np.mod(a+L,N) , np.mod(a+1-b,L) , np.mod(a-1-b,L)]]; 
            dE = 2.0*J*S[a]*(v.sum()) + 2.0*H*S[a];
            if dE <= 0:
                S[a] = -S[a];
                E = E + dE;                
            else:
                p = np.exp(-dE/temp);        
                if np.random.uniform(0,1) <= p:
                    S[a] = -S[a];
                    E = E + dE;
            if t >= 100000:
                m[i,t-100000] = S.sum()/N; # magnetizacao normalizada
        ef[i] = E;
        E = Eini;

        ## Copia condicao inicial para todas realizacoes acontecerem com a mesma condicao inicial ##
        S = Sini.copy();
        ####

    vm = m.std(axis = 0)**2;
    mm = m.mean(axis = 0)**2;
    m.fill(0);
    plt.subplot(221);
    plt.plot(ef);
    ef.fill(0);
    plt.subplot(222);
    plt.plot(Meitempo,mm);
    mm.fill(0);
    plt.subplot(223);
    plt.plot(Meitempo,vm);
    vm.fill(0);
plt.show();
