# -*- coding: utf-8 -*-
import numpy as np
import pylab as p
import copy

random = np.random

class Crossover:
    """ Baseado em
http://en.wikipedia.org/wiki/Edge_recombination_operator and Pyevolve's
Crossovers.G1DListCrossoverEdge method
    """

    def __init__(self, mama, papa):
        self.mae = mama.copy()
        self.pai = papa.copy()
        self.sis = [0]
        self.bro = [0]

    def get_edges(self,ind):
        edg = {}
        ind_list = ind['dna']
        
        for i in range(len(ind_list)):
            a, b = ind_list[i], ind_list[i-1]
            if a not in edg:
                edg[a] = []
            else:
                edg[a].append(b)
            if b not in edg:
                edg[b] = []
            else:
                edg[b].append(a)
        return edg
    
    def merge_edges(self,edge_a, edge_b):
        edges = {}
        for value, near in edge_a.items():
            for adj in near:
                if (value in edge_b) and (adj in edge_b[value]):
                    edges.setdefault(value, []).append(adj)
        return edges
    
    def get_edges_composite(self, mom, dad):
        mom_edges = self.get_edges(mom)
        dad_edges = self.get_edges(dad)
        return (mom_edges, dad_edges, self.merge_edges(mom_edges, dad_edges))

    def cross(self):
        
        mom_edges, dad_edges, merged_edges = self.get_edges_composite(self.mae, self.pai)
        
        for c, u in (self.sis, set(self.mae['dna'])),(self.bro, set(self.pai['dna'])):  
            curr = None            
            tamanho = len(u)
            for i in range(tamanho):                
                curr = random.choice(tuple(u)) if not curr else curr
                c.append(curr)
                u.remove(curr)                
                d = [v for v in merged_edges.get(curr, []) if v in u]
                if d:
                    curr = random.choice(d)
                else:
                    s = [v for v in mom_edges.get(curr, []) if v in u]
                    s += [v for v in dad_edges.get(curr, []) if v in u]
                    curr = random.choice(s) if s else None

    # garante que sempre haverá um 0 no início do dna
        pos0sis = self.sis.index(0)
        self.sis[pos0sis] = self.sis[0]
        self.sis[0] = 0
        pos0bro = self.bro.index(0)
        self.bro[pos0bro] = self.bro[0]
        self.bro[0] = 0
    
        sis0 = self.mae.copy()
        bro0 = self.pai.copy()
        sis0['dna'] = self.sis
        bro0['dna'] = self.bro
    
        return (sis0, bro0)

def mutate(guy):
    """ Mutation by sublist reversing """
    inicio = random.choice(range(1,len(guy['dna'])-1))
    fim = random.choice(range(inicio, len(guy['dna'])-1))
    # invertemos a ordem dessa sublista
    aux = guy.copy()
    foo = aux['dna'][inicio:fim+1]
    foo = foo[::-1]
    # trocamos a sublista antiga pela invertida
    aux['dna'][inicio:fim+1] = foo[:]
    return aux

# PCA #########################################################################

def pca(data):
    data = np.array(data, dtype=float)
    
    # normalizamos a matriz de dados (X = X - mean) e dividimos pelo d.p.
    # X = (X - mean) / dp
    for i in range(data.shape[1]):
        # adiciono um valor irrisorio 0.001 no denominador para nao
        # dar divisao por zero
        data[:,i] = (data[:,i] - data[:,i].mean())/(data[:,i].std()+0.001)
    
    # calculamos a matriz de covariância de X
    matriz_cov = np.cov(data, bias=1, rowvar=0)
    
    # calculamos os autovetores e autovalores e ordenamos em ordem decresc.
    autovalores, autovetores = np.linalg.eig(matriz_cov)
    args = np.argsort(autovalores)[::-1]
    autovalores = autovalores[args]
    autovetores = autovetores[args]
    
    # calculamos os componentes principais para todos os dados
    dados_finais = np.dot(autovetores.T, data.T)
    principais = dados_finais.T
    
    return principais

# Cluster #####################################################################

# Util ########################################################################

def ellipse(x, y, a, b, angle, steps):
    """ Returns the points for an ellipse centered in x, y with size a, b """
    beta = -angle * (np.pi / 180)
    sinbeta = np.sin(beta)
    cosbeta = np.cos(beta)

    alpha = np.linspace(0, 360, steps).T * (np.pi / 180)
    sinalpha = np.sin(alpha)
    cosalpha = np.cos(alpha)

    X = x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta)
    Y = y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta)

    ell = []
    for i in range(steps):
        ell.append([X[i], Y[i]])

    return np.array(ell)

class GeneticAlgorithm:

    def __init__(self, num_cidades = 50, geracoes = 1000,
                 size_pop_ini= 50, taxa_elite = 0.2, optaxa = 0.6, taxa_crossover = 0.6,
                 taxa_muta = 0.3, vizinhos_muta = 0.2, k = 5):
        # Valores iniciais
        self.num_cidades = num_cidades
        self.geracoes = geracoes
        self.size_pop_ini = size_pop_ini
        # number of neighbors in order to do K-mean training and classification of the output
        self.k = k
        # taxas mutacao
        self.taxa_elite = taxa_elite
        self.optaxa = optaxa
        self.taxa_crossover = taxa_crossover
        self.taxa_muta = taxa_muta
        self.vizinhos_muta = vizinhos_muta
        self.population_init()
        self.gera_cidades()

    def new_guy(self):
        dna = np.arange(1, self.num_cidades)
        random.shuffle(dna)        
        d = {'dna': [0] + dna,
             'fitness': .0,
             'score': .0,
             'parents': []}
    
        return d


    def gera_cidades(self):
        self.cidades = ellipse(0, 0, 1, 1, 0, self.num_cidades+1)

#clusters = []

# inicializa a população com novos indivíduos aleatórios

    def population_init(self):
        self.pops = []
        self.pop = []
        for i in range(self.size_pop_ini):
            self.pop.append(self.new_guy())
            
    def fitness(self, guy):
        return 1. / guy['score']

    def dist(self, c1, c2):
        p1 = self.cidades[c1]
        p2 = self.cidades[c2]
        return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

    def score(self, guy):
        """ Scoring based on the sum of distances of all valid paths 
        """
           #return n.abs(p2 - p1)    
        s = 0
        for i in range(len(guy['dna'])-1):
            c1 = guy['dna'][i]
            c2 = guy['dna'][i+1]
            s = s + self.dist(c1, c2)
    
        return s

    def run(self):

        num_elite = int(self.size_pop_ini * self.taxa_elite)
        num_op = int(self.size_pop_ini * self.optaxa)
        num_rand = int(self.size_pop_ini - num_elite - num_op)

        for generation in range(self.geracoes):

           # copia a elite ('num_elite' primeiros) para a nova população
            self.elite = []
            self.new_pop = []
            for i in range(num_elite):
                self.elite.append(copy.deepcopy(self.pop[i]))
                self.new_pop.append(copy.deepcopy(self.pop[i]))

            # aplica operadores de crossover e mutação apenas na elite, criando novos
            for i in range(int(num_op/2)):
                # crossover: aplicado a dois elementos da elite, escolhidos ao acaso                
                mom = random.choice(self.elite)
                dad = random.choice(self.elite)
                sis = None
                bro = None
                if random.uniform() < self.taxa_crossover:
                            
                    (sis, bro) = Crossover(mom, dad).cross()
                else:
                    sis = copy.deepcopy(mom)
                    bro = copy.deepcopy(dad)

                # mutation
                if random.uniform() < self.taxa_muta:
                    sis = mutate(sis)                    
                    bro = mutate(bro)                    

        # store parents
        #sis['parents'] = [dad, mom]
        #bro['parents'] = [mom, dad]
        
                # store new guys in the new pop
                self.new_pop.append(sis)
                self.new_pop.append(bro)

            # o restante de new pop é obtido criando-se novos aleatórios
            for i in range(num_rand):
                ne = self.new_guy()
                self.new_pop.append(ne)

            # calcula o custo de cada indivíduo
            for i in range(self.size_pop_ini):
                sc = self.score(self.new_pop[i])
                self.new_pop[i]['score'] = sc

            # atualiza o fitness de cada indivíduo
            for i in range(self.size_pop_ini):
                fi = self.fitness(self.new_pop[i])
                self.new_pop[i]['fitness'] = fi

            # sobrescreve a população antiga pela nova
            self.pop = self.new_pop[:]
            self.size_pop_ini = len(self.pop)
            
    # *** TODO *** clusteriza
    # components = pca([guy['dna'] for guy in pop])
    # points = []
    # for i in xrange(len(components)):
    # coords = components[i][:2]
    # ref = pop[i]
    # points.append(Point(coords, ref))
    # clusters = kmeans(points, 2, .05)
    # escolhe representante
     
    # ordenamos a população em função de seu fitness (maior -> menor)
            self.pop.sort(key=lambda x: x['fitness'], reverse=True)
            self.pops.append({'generation': generation,
                              'pop': self.pop,
                              'best': self.pop[0],
                              'best fitness': self.pop[0]['fitness'],
                              'fitness avg': sum([x['fitness'] for x in self.pop]) / len(self.pop),
                              'fitness min': min([x['fitness'] for x in self.pop]),
                              'fitness max': max([x['fitness'] for x in self.pop]),
                              'score avg': sum([x['score'] for x in self.pop]) / len(self.pop),
                              'score min': min([x['score'] for x in self.pop]),
                              'score max': max([x['score'] for x in self.pop])})
          

            print ('*' * 10)
            print('generation: ',generation)
            print('best ',self.pop[0]['dna'])
            print('best fitness: ',self.pop[0]['fitness'])
            print('max fitness: ',max([x['fitness'] for x in self.pop])) 



            
model = GeneticAlgorithm(num_cidades = 40, geracoes = 1000,
                         size_pop_ini= 50, taxa_elite = 0.2, optaxa = 0.6, taxa_crossover = 0.6,
                         taxa_muta = 0.3, vizinhos_muta = 0.2, k = 5)
model.run()


###### PLOT MODEL ######

p.figure()
x = []
y = []
yerr_max = []
yerr_min = []
ymin = []
ymax = []
for po in model.pops:
    x.append(po['generation'])
    y.append(po['fitness avg'])
    ymin.append(po['fitness min'])
    ymax.append(po['fitness max'])
    yerr_max.append(ymax)
    yerr_min.append(ymin)

p.plot(x, ymin, 'g-')
p.plot(x, ymax, 'r-')
p.plot(x, y, '-')
p.xlabel('Generation')
p.ylabel('Fitness Min/Avg/Max')
p.grid(True)

p.figure()
x = []
y = []
yerr_max = []
yerr_min = []
ymin = []
ymax = []
for po in model.pops:
    x.append(po['generation'])
    y.append(po['score avg'])
    ymin.append(po['score min'])
    ymax.append(po['score max'])
    yerr_max.append(ymax)
    yerr_min.append(ymin)

p.plot(x, ymin, 'g-')
p.plot(x, ymax, 'r-')
p.plot(x, y, '-')
p.xlabel('Generation')
p.ylabel('Score Min/Avg/Max')
p.grid(True)


p.show()
