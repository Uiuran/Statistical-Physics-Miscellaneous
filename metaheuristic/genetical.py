# -*- coding: utf-8 -*-

import numpy as np
import random as r
import pylab as p
import copy

# Valores iniciais
cidades = 50
geracoes = 1000
tamanho_populacao = 200

# numero de votantes para decisao
knum = 3

# taxas mutacao
taxa_elite = 0.2
optaxa = 0.6
taxa_crossover = 0.6
taxa_muta = 0.1
vizinhos_muta = 0.1

def crossover(mae, pai):
    """ Baseado em
http://en.wikipedia.org/wiki/Edge_recombination_operator and Pyevolve's
Crossovers.G1DListCrossoverEdge method
    """
    gmae = mae.copy()
    gpai = pai.copy()
    irma = []
    irmao = []

    def get_edges(ind):
        edg = {}
        ind_list = ind['dna']
        
        for i in xrange(len(ind_list)):
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
    
    def merge_edges(edge_a, edge_b):
        edges = {}
        for value, near in edge_a.items():
            for adj in near:
                if (value in edge_b) and (adj in edge_b[value]):
                    edges.setdefault(value, []).append(adj)
        return edges
    
    def get_edges_composite(mom, dad):
        mom_edges = get_edges(mom)
        dad_edges = get_edges(dad)
        return (mom_edges, dad_edges, merge_edges(mom_edges, dad_edges))
    
    mom_edges, dad_edges, merged_edges = get_edges_composite(gmae, gpai)

    for c, u in (irma, set(gmae['dna'])), (irmao, set(gpai['dna'])):
        curr = None
        for i in xrange(len(gmae['dna'])):
            curr = r.choice(tuple(u)) if not curr else curr
            c.append(curr)
            u.remove(curr)
            d = [v for v in merged_edges.get(curr, []) if v in u]
            if d:
                curr = r.choice(d)
            else:
                s = [v for v in mom_edges.get(curr, []) if v in u]
                s += [v for v in dad_edges.get(curr, []) if v in u]
                curr = r.choice(s) if s else None
    # garante que sempre haverá um 0 no início do dna
    pos0irma = irma.index(0)
    irma[pos0sister] = irma[0]
    irma[0] = 0
    pos0brother = irmao.index(0)
    irmao[pos0brother] = irmao[0]
    irmao[0] = 0
    
    irma0 = gmae.copy()
    irmao0 = gpai.copy()
    irma0['dna'] = irma
    irmao0['dna'] = irmao
    
    return (irma0, irmao0)

def mutate(guy):
    """ Mutation by sublist reversing """
    inicio = r.choice(range(1,len(guy['dna'])-1))
    fim = r.choice(range(inicio, len(guy['dna'])-1))
    # invertemos a ordem dessa sublista
    aux = guy.copy()
    foo = aux['dna'][inicio:fim+1]
    foo.reverse()
    # trocamos a sublista antiga pela invertida
    aux['dna'][inicio:fim+1] = foo[:]
    return aux

def fitness(guy):
    return 1. / guy['score']

def score(guy):
    """ Scoring based on the sum of distances of all valid paths """
    def dist(c1, c2):
        p1 = cidades[c1]
        p2 = cidades[c2]
        return n.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
        #return n.abs(p2 - p1)
    
    s = 0
    for i in xrange(len(guy['dna'])-1):
        c1 = guy['dna'][i]
        c2 = guy['dna'][i+1]
        s = s + dist(c1, c2)
    
    return s

def new_guy():
    dna = range(1, cidades)
    r.shuffle(dna)
    
    d = {'dna': [0] + dna,
         'fitness': .0,
         'score': .0,
         'parents': []}
    return d.copy()

# PCA #########################################################################

def pca(data):
    data = np.array(data, dtype=float)
    
    # normalizamos a matriz de dados (X = X - mean) e dividimos pelo d.p.
    # X = (X - mean) / dp
    for i in xxrange(data.shape[1]):
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
    for i in xrange(steps):
        ell.append([X[i], Y[i]])

    return np.array(ell)

# Main loop ###################################################################

num_elite = int(tamanho_populacao * taxa_elite)
num_op = int(tamanho_populacao * optaxa)
num_rand = int(tamanho_populacao - num_elite - num_op)
print num_elite, num_op, num_rand
cidades = ellipse(0, 0, 1, 1, 0, cidades+1)

clusters = []

# inicializa a população com novos indivíduos aleatórios
pops = []
pop = []
for i in xrange(tamanho_populacao):
    pop.append(new_guy())

for generation in xrange(geracoes):
    # copia a elite ('num_elite' primeiros) para a nova população
    elite = []
    new_pop = []

    for i in xrange(num_elite):
        elite.append(copy.deepcopy(pop[i]))
        new_pop.append(copy.deepcopy(pop[i]))

    # aplica operadores de crossover e mutação apenas na elite, criando novos
    for i in xrange(num_op/2):
        # crossover: aplicado a dois elementos da elite, escolhidos ao acaso
        mom = r.choice(elite)
        dad = r.choice(elite)
        sis = None
        bro = None
        if r.random() < taxa_crossover:
            (sis, bro) = crossover(mom, dad)
        else:
            sis = copy.deepcopy(mom)
            bro = copy.deepcopy(dad)

        # mutation
        if r.random() < taxa_muta:
            sis = mutate(sis)
            bro = mutate(bro)

        # store parents
        #sis['parents'] = [dad, mom]
        #bro['parents'] = [mom, dad]
        
        # store new guys in the new pop
        new_pop.append(sis)
        new_pop.append(bro)

    # o restante de new pop é obtido criando-se novos aleatórios
    for i in xrange(num_rand):
        ne = new_guy()
        new_pop.append(ne)

    # calcula o custo de cada indivíduo
    for i in xrange(tamanho_populacao):
        sc = score(new_pop[i])
        new_pop[i]['score'] = sc

    # atualiza o fitness de cada indivíduo
    for i in xrange(tamanho_populacao):
        fi = fitness(new_pop[i])
        new_pop[i]['fitness'] = fi

    # sobrescreve a população antiga pela nova
    pop = new_pop[:]

    # clusteriza
    # components = pca([guy['dna'] for guy in pop])
    # points = []
    # for i in xrange(len(components)):
    # coords = components[i][:2]
    # ref = pop[i]
    # points.append(Point(coords, ref))
    # clusters = kmeans(points, 2, .05)

    # escolhe representante

    # *** TODO ***

    # ordenamos a população em função de seu fitness (maior -> menor)
    pop.sort(key=lambda x: x['fitness'], reverse=True)

    pops.append({'generation': generation,
                 'pop': pop,
                 'best': pop[0],
                 'best fitness': pop[0]['fitness'],
                 'fitness avg': sum([x['fitness'] for x in pop]) / len(pop),
                 'fitness min': min([x['fitness'] for x in pop]),
                 'fitness max': max([x['fitness'] for x in pop]),
                 'score avg': sum([x['score'] for x in pop]) / len(pop),
                 'score min': min([x['score'] for x in pop]),
                 'score max': max([x['score'] for x in pop])})


    print '*' * 10
    print 'generation: ', generation
    print 'best ', pop[0]['dna']
    print 'best fitness: ', pop[0]['fitness']
    print 'max fitness: ', max([x['fitness'] for x in pop])


# Fitness

p.figure()
x = []
y = []
yerr_max = []
yerr_min = []
ymin = []
ymax = []
for po in pops:
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
for po in pops:
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
