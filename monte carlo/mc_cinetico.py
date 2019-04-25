import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import Image as img

T = 1000; #numero total de passos do algoritmo
n0 = 10; #numero inicial de individuos

def mc_cinetico(pop0 , Tt, taxa_nascimento0, taxa_morte0):
    ''' Algoritmo Monte Carlo Cinetico, exemplificado atraves do processo de nascimento e morte '''
    taxa_global = 1.0;
    n = pop0; # numero de individuos
    taxa_nascimento = taxa_nascimento0;
    taxa_morte = taxa_morte0;
    dT = 0.0;
    time = [];
    pop = [];
    for i in range(Tt):
        time.append(dT);
        pop.append(n);
        if n == 0:
            break;
  
        taxa_global = n*(taxa_nascimento + taxa_morte);                  
        dT = dT + np.random.exponential(1.0/taxa_global);
        eps = np.random.uniform(0.0,1.0);
        frac_nascimento = taxa_nascimento/taxa_global;
        frac_morte = taxa_morte/taxa_global;
        if eps < n*frac_nascimento:
            n = n + 1;
        else:
            n = n - 1;
    return time,pop;

def taxa_contato(populacao, taxa_autocatalise, L):
    ''' Calcula taxa global de probabilidades para o processo de contato (com 2 estados individuais: 1 ou 0) em rede regular 2D. '''
    N = len(populacao);
    taxas_individuais = np.zeros(N);
    for a in range(N):
        b = np.floor(a/L)*L;
        v = populacao[[np.mod(a-L,N) , np.mod(a+L,N) , b+np.mod(a+1-b,L) , b+np.mod(a-1-b,L) ]];
        taxas_individuais[a] = 0.25*taxa_autocatalise*(1.0-populacao[a])*v.sum() + populacao[a];
    taxa_global = taxas_individuais.sum();
    return taxa_global, taxas_individuais;

def mc_contato(populacao, Tt, taxa_autocatalise, L, out = 'arquivo', experiment_title = '0'):
    ''' Algoritmo Monte Carlo Cinetico para processo de contato. Opcoes de saidas (out) sao: 'arquivo' grava simulacao, sequencia de {0,1}, em arquivo .ssd (stochastic simulation data), 'filme' para gravar simulacao como uma serie de frames da grade 2D e 'tela' para imprimir a simulacao no prompt de comando'''
    
    N = len(populacao); # numero de individuos  
    dT = 0.0;
    tstep = 0; # contador de passos temporais
    time = [];
    taxa_global, taxas_individuais = taxa_contato(populacao, taxa_autocatalise, L); #calcula taxa global e as individuais do estado inicial

    beta = lambda x,eps: np.floor(x/eps)*eps;
    alpha = lambda x,eps,delta,b: [np.mod(x -eps,delta) , np.mod(x +eps,delta) , b+np.mod(x+1-b,eps) , b+np.mod(x-1-b,eps)];
    if Tt > 10.0:
        Tt = 10.0;   
   
    while Tt > 0.0:
                
        if out == 'filme':
            d = [];
            for y in range(L):
                for j in range(L):
                    b0 = y*L+j
#                    b = beta(b0,L);
#                    l = populacao[alpha(b0,L,N,b)];
#                    l = l.tolist();
                    if populacao[b0] == 0:
                        d.append((000,255,000));
                    else:
                        d.append((255,000,000));
            pic = img.new('RGB',(L,L));
            pic.putdata(d);        
            pic.save('frames/c'+str(tstep+1)+'.png');                 
       
        probabilidade_cumulativa = taxas_individuais.cumsum()/taxa_global;               
        a = np.random.exponential(1.0/taxa_global);
        Tt = Tt - a;
        dT = dT +  a;       
        transicao = np.random.uniform(0.0,1.0);
        loctransicao = probabilidade_cumulativa.searchsorted(transicao);
        
        tstep = tstep + 1;
        time.append(dT);

        if loctransicao == L:
            loctransicao = loctransicao - 1;

        if populacao[loctransicao] == 1:
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);

            taxa_global = taxa_global - (taxas_individuais[loctransicao] + taxas_individuais[l].sum());
            populacao[loctransicao] = 0;                 
            taxas_individuais[loctransicao] = 0.25*taxa_autocatalise*populacao[l].sum();

	    a = np.mod(loctransicao-L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];

            a = np.mod(loctransicao+L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];

            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao+1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];

            
            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao-1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);            
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];
                
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);            
         
            taxa_global = taxa_global + taxas_individuais[loctransicao] + taxas_individuais[l].sum();
        else:
            
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);

            taxa_global = taxa_global - (taxas_individuais[loctransicao] + taxas_individuais[l].sum());

            populacao[loctransicao] = 1;                 
            taxas_individuais[loctransicao] = populacao[loctransicao];
             
            a = np.mod(loctransicao-L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);             
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];

            a = np.mod(loctransicao+L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];

            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao+1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];

            
            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao-1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            taxas_individuais[a] = 0.25*taxa_autocatalise*(1-populacao[a])*populacao[l].sum() + populacao[a];
                
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);
            taxa_global = taxa_global + taxas_individuais[loctransicao] + taxas_individuais[l].sum();
            
        if out == 'tela':
            print '-----------------';
            for y in range(L):           
                print populacao[y*L:(y+1)*L];
            print '*****************';
        if out == 'arquivo':
            np.savetxt('data/'+experiment_title+'_'+str(tstep)+'.ssd', populacao, fmt='%-10.1i');
    np.savetxt('data/'+experiment_title+'_time.t', time, fmt='%-10.8f');

def movie(exptitle, time_num, L, tipo='contato'):

    beta = lambda x,eps: np.floor(x/eps)*eps;
    alpha = lambda x,eps,delta,b: [np.mod(x -eps,delta) , np.mod(x +eps,delta) , b+np.mod(x+1-b,eps) , b+np.mod(x-1-b,eps)];
    if tipo == 'contato':
        for i in range(time_num):
            d = [];
            populacao = np.loadtxt('./data/'+exptitle+'_'+str(i+1)+'.ssd');
            for y in range(L):
                for j in range(L):
                    b0 = y*L+j
#                b = beta(b0,L); ##caso seja necessario informacao de primeiros vizinhos de b0 na rede regular (dado por beta e alpha)
#                l = populacao[alpha(b0, L, L*L,b)];
#                l = l.tolist();
                    if populacao[b0] == 0:
                        d.append((000,255,000));
                    else:
                        d.append((255,000,000));
            pic = img.new('RGB',(L,L));
            pic.putdata(d);        
            pic.save('frames/c'+str(i)+'.png');
    elif tipo == 'predador presa':
        tstep = 0;
        for i in range(time_num):
            d = [];
            populacao = np.loadtxt('./data/'+exptitle+'_'+str(i+1)+'.ssd');
            for y in range(L):
                for j in range(L):
                    b0 = y*L+j
#                    b = beta(b0,L);
#                    l = populacao[alpha(b0,L,N,b)];
#                    l = l.tolist();
                    if populacao[b0] == 0:
                        d.append((000,255,000));
                    elif populacao[b0] == 1:
                        d.append((255,000,000));
                    else:
                        d.append((255,255,000));

            pic = img.new('RGB',(L,L));
            pic.putdata(d);        
            pic.save('frames/pp'+str(tstep+1)+'.png');
            tstep = tstep + 1;

   

def taxa_predadorpresa(populacao, taxas, L):
    ''' Calcula taxa global de probabilidades para o processo de contato (com 2 estados individuais: 1 ou 0) em rede regular 2D. '''
    N = len(populacao);
    taxas_individuais = np.zeros(N);
    if len(taxas) != 3:
        print 'devem haver 3 taxas para este sistema';
        return 0;

    for a in range(N):
        b = np.floor(a/L)*L;
        v = populacao[[np.mod(a-L,N) , np.mod(a+L,N) , b+np.mod(a+1-b,L) , b+np.mod(a-1-b,L) ]];
        vv = v.tolist();
        taxas_individuais[a] = 0.25*(taxas[0]*[populacao[a]].count(0)*vv.count(1) + taxas[1]*[populacao[a]].count(1)*vv.count(2)) + taxas[2]*[populacao[a]].count(2);
    taxa_global = taxas_individuais.sum();
    return taxa_global, taxas_individuais;



def mc_predadorpresa(populacao, Tt, taxa_autocatalise, L, out = 'arquivo', experiment_title = '0'):
    ''' Algoritmo Monte Carlo Cinetico para processo tipo Predador-Presa. Opcoes de saidas (out) sao: 'arquivo' grava simulacao, sequencia com numeros {0,1,2}, em arquivo .ssd (stochastic simulation data), 'filme' para gravar simulacao como uma serie de frames da grade 2D e 'tela' para imprimir a simulacao no prompt de comando'''
    
    N = len(populacao); # numero de individuos  
    dT = 0.0;
    tstep = 0; # contador de passos temporais
    time = [];
    taxa_global, taxas_individuais = taxa_predadorpresa(populacao, taxa_autocatalise, L); #calcula taxa global e as individuais do estado inicial

    beta = lambda x,eps: np.floor(x/eps)*eps;
    alpha = lambda x,eps,delta,b: [np.mod(x -eps,delta) , np.mod(x +eps,delta) , b+np.mod(x+1-b,eps) , b+np.mod(x-1-b,eps)];
    if Tt > 10.0:
        Tt = 10.0;   
   
    while Tt > 0.0:
                
        if out == 'filme':
            d = [];
            for y in range(L):
                for j in range(L):
                    b0 = y*L+j
#                    b = beta(b0,L);
#                    l = populacao[alpha(b0,L,N,b)];
#                    l = l.tolist();
                    if populacao[b0] == 0:
                        d.append((000,255,000));
                    elif populacao[b0] == 1:
                        d.append((255,000,000));
                    else:
                        d.append((255,255,000));

            pic = img.new('RGB',(L,L));
            pic.putdata(d);        
            pic.save('frames/'+str(tstep+1)+'.png');                 
       
        probabilidade_cumulativa = taxas_individuais.cumsum()/taxa_global;               
        a = np.random.exponential(1.0/taxa_global);
        Tt = Tt - a;
        dT = dT +  a;       
        transicao = np.random.uniform(0.0,1.0);
        loctransicao = probabilidade_cumulativa.searchsorted(transicao);
        
        tstep = tstep + 1;
        time.append(dT);

        if loctransicao == L:
            loctransicao = loctransicao - 1;

        if populacao[loctransicao] == 1:
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);

            taxa_global = taxa_global - (taxas_individuais[loctransicao] + taxas_individuais[l].sum());
            populacao[loctransicao] = 2;                 
            taxas_individuais[loctransicao] = taxa_autocatalise[2];

	    a = np.mod(loctransicao-L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);

            a = np.mod(loctransicao+L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);

            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao+1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);

            
            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao-1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);
                
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);            
         
            taxa_global = taxa_global + taxas_individuais[loctransicao] + taxas_individuais[l].sum();
        elif populacao[loctransicao] == 0:
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);
            v = populacao[l];
            vv = v.tolist();

            taxa_global = taxa_global - (taxas_individuais[loctransicao] + taxas_individuais[l].sum());
            populacao[loctransicao] = 1;                 
            taxas_individuais[loctransicao] = 0.25*taxa_autocatalise[1]*vv.count(2);
             
            a = np.mod(loctransicao-L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);             
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);
 

            a = np.mod(loctransicao+L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);


            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao+1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);
            
            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao-1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);
                
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);
            taxa_global = taxa_global + taxas_individuais[loctransicao] + taxas_individuais[l].sum();
        else:
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);

            taxa_global = taxa_global - (taxas_individuais[loctransicao] + taxas_individuais[l].sum());
            populacao[loctransicao] = 0;  
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[loctransicao] = 0.25*taxa_autocatalise[0]*vv.count(1);
            
            a = np.mod(loctransicao-L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);

            a = np.mod(loctransicao+L,N);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);

            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao+1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);
            
            b = beta(loctransicao,L);
            a = b+np.mod(loctransicao-1-b,L);
            b = beta(a,L);
            l = alpha(a,L,N,b);
            v = populacao[l];
            vv = v.tolist();
            taxas_individuais[a] = 0.25*(taxa_autocatalise[0]*[populacao[a]].count(0)*vv.count(1) + taxa_autocatalise[1]*[populacao[a]].count(1)*vv.count(2)) + taxa_autocatalise[2]*[populacao[a]].count(2);
                
            b = beta(loctransicao,L);
            l = alpha(loctransicao,L,N,b);
            taxa_global = taxa_global + taxas_individuais[loctransicao] + taxas_individuais[l].sum();       
        if out == 'tela':
            print '-----------------';
            for y in range(L):           
                print populacao[y*L:(y+1)*L];
            print '*****************';
        if out == 'arquivo':
            np.savetxt('data/'+experiment_title+'_'+str(tstep)+'.ssd', populacao, fmt='%-10.1i');
    np.savetxt('data/'+experiment_title+'_time.t', time, fmt='%-10.8f');

def encodaMovie(diretorio = 'frames/', nome = 'saida'):
    ''' Faz filme de arquivos imagens a partir de 'diretorio', com nome 'nome' ''' 
    string_linha_de_commando = 'cd '+diretorio+' && mencoder -noskip mf://*.png -mf fps=20:type=png -ovc lavc -lavcopts vcodec=mpeg4 -o '+nome+'.avi'
    os.system(string_linha_de_commando);



#####################################################################################################
#pop = np.random.uniform(0,1, 10000);
#pop0 = np.array([1 if a < 0.5 else 0 for a in pop]);
#popp = np.random.randint(0,3,10000);
#mc_predadorpresa(popp,1.5, [1.756,2.85,0.75], 100, 'arquivo', experiment_title = 'pp1');
#movie('pp1', 9105 , 100, tipo = 'predador presa');
encodaMovie();
