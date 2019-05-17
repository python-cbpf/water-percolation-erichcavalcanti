def test(ord):
    """ Só plota"""
    Px,Py = Clust_see(21,10,100,True)
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(Px,Py[ord],marker='o',ls=':')
    ax.set_xlabel('Free cell probability')
#    ax.set_ylabel('Average number of clusters')
#    print(Py)
    return

def test2(L1,L2,L3,cut):
    """ Só plota"""
    Px1,Py1 = Clust_see(L1,10,100,cut)
    Px2,Py2 = Clust_see(L2,10,100,cut)
    Px3,Py3 = Clust_see(L3,10,100,cut)

    ax1 = plt.subplot(311)
    plt.plot(Px1,Py1[0],'ro-')
    plt.plot(Px2,Py2[0],'go-')
    plt.plot(Px3,Py3[0],'bo-')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax2 = plt.subplot(312, sharex=ax1)
    plt.plot(Px1,Py1[1],'ro-')
    plt.plot(Px2,Py2[1],'go-')
    plt.plot(Px3,Py3[1],'bo-')
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax3 = plt.subplot(313, sharex=ax1)
    plt.plot(Px1,Py1[2],'ro-')
    plt.plot(Px2,Py2[2],'go-')
    plt.plot(Px3,Py3[2],'bo-')
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.show()
    
    return

def Clust_see(L,prob_points,samples,tocut):
    """ Faz Clust_Avg para um conjunto de probabilidades pf. 
    Organiza dados para facilitar plot"""
    pf_range = np.linspace(0,1,prob_points)
    l_nl = np.array([Clust_Avg(L,pf,samples,tocut) for pf in pf_range])
    return pf_range, l_nl.T

from collections import Counter
def Clust_Avg(L,pf,samples,tocut) :
    """ Após obter NS,SNS e S2NS, tira a média para certo conjunto de samples    
    lista de todos os clusters (com seus tamanhos)
        CLUST = Get_Clusters(L,pf)
        CLUST[0]
    lista de todos os clusters que não tocam a borda (e seus tamanhos)
        CLUST_CUT = np.array(CLUST[0])[np.nonzero(CLUST[2])]
    reorganização da lista, por contagem
        CLUST_LIST_CUT = Counter(CLUST_CUT)
        CLUST_LIST = Counter(CLUST)
    
    N_S numero de clusters com tamanho S dentro do grid de tamanho L
    N_S_CUT removendo da contagem os que tocam a borda
    \sum_S N_S = número total de clusters presente
    
    SNS = S * N_S, S2NS = S^2 * N_S    
    """
    if tocut : NS_CUT = 0; SNS_CUT = 0; S2NS_CUT = 0
    else : NS = 0; SNS = 0; S2NS = 0    

    for i in range(samples) : 
        """ Trata as 'estatisticas' do Get_Clusters obtendo NS,SNS e S2NS
        para o caso completo e para o caso onde os clusters que tocam na borda
        são ignorados (CUT)"""
        CLUST = Get_Clusters(L,pf)
 
        if tocut :
            CLUST_CUT = np.array(CLUST[0])[np.nonzero(CLUST[2])] #remove da lista os clusters que tocam borda
            NS_CUT += len(CLUST_CUT) #número total de clusters        
            CLUST_LIST_CUT = Counter(CLUST_CUT)
            if len(list(enumerate(CLUST_LIST_CUT)))>0 :
                SNS_CUT_aux, S2NS_CUT_aux = sum( np.array([ S[1]*CLUST_LIST_CUT[S[1]], S[1]*S[1]*CLUST_LIST_CUT[S[1]] ]) for S in enumerate(CLUST_LIST_CUT))
                SNS_CUT += SNS_CUT_aux; S2NS_CUT += S2NS_CUT_aux 
        else :
            NS += len(CLUST[0]) #número total de clusters
            CLUST_LIST = Counter(CLUST[0])
            if len(list(enumerate(CLUST_LIST)))>0 :
                SNS_aux, S2NS_aux = sum( np.array([ S[1]*CLUST_LIST[S[1]], S[1]*S[1]*CLUST_LIST[S[1]] ]) for S in enumerate(CLUST_LIST))
                SNS += SNS_aux; S2NS += S2NS_aux

    if tocut :
        NS_CUT /= samples; SNS_CUT /= samples;  S2NS_CUT /= samples 
        return NS_CUT,SNS_CUT,S2NS_CUT
    else :
        NS /= samples; SNS /= samples; S2NS /= samples;
        return NS,SNS,S2NS
    return

def Get_Clusters(G_L,G_p):
    """ Dado um tamanho de rede e probabilidade de célula livre, retorna dados
    dos clusters presentes: número de conexões de cada cluster, grid inicial,
    informação se os clusters tocam as bordas do grid"""
    ##inicialização
#    G_L = 31
#    G_p = .6
    x = lattice(G_L,G_p)
    Fx = x.FREE_GRID
    ##conjunto de clusters
    Clust = []
    Clust_conn = []
    Clust_edges = []
    ##Primeiro Cluster
    x.percolatefull()
    Clust.append(x.GRID) # armazena cluster
    Clust_conn.append(Clust[len(Clust_conn)].sum()) #numero de conexões do cluster
#        #conferia se tocava a borda pelo menos uma vez
#    x.checkedge()
#    Clust_edges.append(not(x.edgeB))
    #confere se cluster toca a si mesmo pelas bordas
    Clust_edges.append(not(x.checkedgeperc())) 
    FlagCluster = x.hasCluster() # confere se ainda existe algum cluster a ser encontrado
    
    while FlagCluster:
        x.restart(Clust[len(Clust_conn)-1])    
        x.percolatefull()
        Clust.append(x.GRID) # armazena cluster
        Clust_conn.append(Clust[len(Clust_conn)].sum()) #numero de conexões do cluster
#        #conferia se tocava a borda pelo menos uma vez
#        x.checkedge()
#        Clust_edges.append(not(x.edgeB)) #fará filtro, só sobrevive quem não toca borda
        #confere se cluster toca a si mesmo pelas bordas
        Clust_edges.append(not(x.checkedgeperc())) #fará filtro, só sobrevive quem não se toca

        FlagCluster = x.hasCluster() # confere se ainda existe algum cluster a ser encontrado
      
#    ##Conferindo se obtivemos realmente todos os clusters
#    All_Clust = Clust[0].copy()
#    for i in range(1,len(Clust_conn)):
#        All_Clust += Clust[i]
#    0 in (Fx == All_Clust)
    return Clust_conn,Fx,Clust_edges

###############################################################################

def just_do_it():
    
    L_range = ([51,101,151,201,251])
#    L_range = ([301,401,501,601,701])
    T_values = []
    
    for L in L_range:
        T_values.append(runover_main(L,1000,10,'A',True,False)) 

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(L_range,T_values,marker='o',ls=':')
    ax.set_xlabel('Grid size (L)')
    ax.set_ylabel('Run time (s)')

    print(L_range)
    print(T_values)
    return

###############################################################################

import numpy as np
from matplotlib import pyplot as plt
import time
from pylab import matplotlib as mpl
from scipy.optimize import curve_fit
mpl.rcParams.update({'font.size': 18, 'font.family': 'serif'})


### -------------- Class : lattice ---------------
class lattice:
    """
    REQUIRES : 
        import numpy as np
        from matplotlib import pyplot as plt
    INPUT
        Grid size : L0 -> self.L
        Free cell probability : pf
    MIDDLE OUTPUT
        FREE_GRID : path to follow
    OUTPUT
        GRID / NEW_GRID; edge / edgeB, edgeU, edgeD, edgeL, edgeR; conv; STEP; density(); ? fractaldimension : outside function that uses self.toplot()
    DESCRIPTION
        * We use a GRID with size (L+2)^2. The border of the GRID is neglected.
            Input and Output does not seem the "extra" size
        * Comment: 
            It is assumed that the border is fully free. 
            at some moment we can explore some boundary condition
        * To simplify we use odd L, so that the grid center is well defined
    """   
    def __init__(self,L0,pf):
        if (L0%2) : self.L = L0
        else : self.L = L0+1; print("Size changed to "+str(L0+1))
        
        self.GRID = np.zeros((self.L+2,self.L+2),dtype=bool)    
        Lm=np.int((self.L+1)/2)        
        self.GRID[Lm,Lm]=True
        
        INIT_RAND = np.random.rand(self.L,self.L)
        self.FREE_GRID = np.zeros((self.L+2,self.L+2),dtype=bool)    
        self.FREE_GRID[1:-1,1:-1] = (INIT_RAND < pf) #True se vazio
        self.FREE_GRID[Lm,Lm]=True

        self.edgeU=False; self.edgeD=False; self.edgeL=False; self.edgeR=False
        self.edge=False; self.conv=False; self.STEP=0
        return
    def hasCluster(self):
        return (1 in self.FREE_GRID ^ self.GRID)
    def restart(self,CLUSTER):
        """ Útil para obter vários clusters"""
        if self.hasCluster():
            self.FREE_GRID = self.FREE_GRID ^ CLUSTER
    
            # Melhorar isso... .nonzero() dá lista de TODOS! Quero somente o primeiro
            coord = (np.array(self.FREE_GRID.nonzero()).T[0,0],
                     np.array(self.FREE_GRID.nonzero()).T[0,1])
            self.GRID = np.zeros((self.L+2,self.L+2),dtype=bool)    
            self.GRID[coord]=True
               
            self.edgeU=False; self.edgeD=False; self.edgeL=False; self.edgeR=False
            self.edge=False; self.conv=False; self.STEP=0
        else : print("Não existem mais clusters!")
        return        
    def toprint(self):
        return self.GRID[1:-1,1:-1]*1 + np.invert(self.FREE_GRID[1:-1,1:-1])*(-1)
    def nextstep(self):
        """ Método sem condição de contorno periodica. Bordas são ignoradas"""
        self.NEW_GRID = self.GRID.copy()
        hasOne = (self.GRID[0:-2,1:-1] + self.GRID[1:-1,0:-2] + 
                  self.GRID[1:-1,2:] + self.GRID[2:  ,1:-1])
        self.NEW_GRID[1:-1,1:-1] = self.GRID[1:-1,1:-1] + self.FREE_GRID[1:-1,1:-1]*hasOne
        return
    def acept(self):
        self.GRID = self.NEW_GRID
        return
    def checkedgeperc(self):
        """ GRID with size self.L+2
            index 0 and L+1 are null borders
            index 1 and L are true edges
            VERIFY if cluster itself through borders.
        """
        Horizontal = (1 in self.GRID.T[1] & self.GRID.T[self.L])
        Vertical = (1 in self.GRID[1] & self.GRID[self.L])
        return Horizontal | Vertical
    def checkedge(self):
        """ GRID with size self.L+2
            index 0 and L+1 are null borders
            index 1 and L are true edges
        """
        self.edgeU = 1 in self.GRID[1]
        self.edgeD = 1 in self.GRID[self.L]
        self.edgeL = 1 in self.GRID.T[1]
        self.edgeR = 1 in self.GRID.T[self.L]
        self.edge = self.edgeU + self.edgeD + self.edgeL + self.edgeR
        self.edgeB = np.bool(self.edge)
        return 
    def checkconvergence(self):
        self.conv = False not in (self.GRID == self.NEW_GRID)
        return
    def percolate(self):
        MAX_STEPS = 1000
        while (self.STEP < MAX_STEPS):
            self.nextstep()
            self.checkedge()
            self.checkconvergence()
            if self.conv : break
            if self.edge : break
            self.acept()
            self.STEP +=1
        return
    def percolatefull(self):
        """ Removido condição de parada ao chegar na borda 
        mas NÃO apresenta condição de contorno periodica
        """
        MAX_STEPS = 1000
        while (self.STEP < MAX_STEPS):
            self.nextstep()
            self.checkconvergence()
            if self.conv : break
            self.acept()
            self.STEP +=1
        return
    def BC_periodic(self):
        """ Adiciona condições de contorno periodicas"""
        self.GRID[0] = self.GRID[self.L]
        self.GRID[self.L+1] = self.GRID[1]
        self.GRID.T[0] = self.GRID.T[self.L]
        self.GRID.T[self.L+1] = self.GRID.T[1]        
        return
    def percolatefull_periodic(self):
        """ Removido condição de parada ao chegar na borda.
        Adiciona condição de contorno periodica
        """
        MAX_STEPS = 1000
        while (self.STEP < MAX_STEPS):
            self.BC_periodic()
            self.nextstep()
            self.checkconvergence()
            if self.conv : break
            self.acept()
            self.STEP +=1
        return
    def evolve(self,N):
        for i in range(N):
            self.nextstep()
            self.acept()
            self.STEP +=1
        return
    def density(self):
        return self.GRID.sum() / (self.L*self.L)
    def density2(self):
        return self.GRID.sum() / self.FREE_GRID.sum()

### -------------- Function to measure time. Useful for optimization ---------------
import cProfile
def measure(__func__):
    """ medir tempo de execução dos processos internos de uma função """
    if __name__ == '__main__':
        pr = cProfile.Profile()
        pr.enable()
        __func__()
        pr.disable()
        pr.print_stats()
    return

### -------------- Just plot ---------------
def plot(toprint):
    plt.figure()
    plt.imshow(toprint,vmin=-1,vmax=1,cmap='coolwarm')
    plt.axis('off')
    return

################# MODIFICAR : tornar funções mais limpas #################
### -------------- Function : average of some property ---------------
def average_main(L,pf,samples,kind):
    """A diferença do kind é na escolha do Return na função main()
        Isso significa, aqui, a escolha do que estamos tomando a média no número
        de samples. """
        
    if ((kind == 'A') or (kind=='B') or (kind=='E')): average=0
    elif ((kind == 'C') or (kind=='D')): average=np.zeros(4)
    
    for it in range(samples):
        x = lattice(L,pf)
        x.percolate()
        if (kind=='A'): average += x.edgeB
        elif (kind=='B'): average += x.density()
        elif (kind=='C'): average += [x.edgeU,x.edgeD,x.edgeL,x.edgeR]
        elif ((kind=='D') and (x.edge>0)): average[x.edge-1] += 1
        elif (kind=='E') : average += fractal_dimension(x.toprint())
    average /= samples
    
    return average

### -------------- Function : salva dados, faz plot, faz fit ---------------

def runover_main(L,samples,prob_points,kind,FlagFit,FlagPlot):
    """ average_main ao longo de um intervalo de valores de probabilidade
    """
    TheMethod = ['Grid Copies', 'Border Expansion','Grid Copies with Bits']
    method=0
    if (kind=='A'): thestrings = ['Percolation','Frequency']
    elif (kind=='B'): thestrings = ['Density','Mean Density']
    elif (kind=='C'): 
        thestrings = ['PercWhich','Frequency']
        thelabels = ['up','down','left','right']
    elif (kind=='D'): 
        thestrings = ['PercMany','Frequency']
        thelabels = ['1 wall','2 walls','3 walls','4 walls']
    elif (kind=='E'): 
        thestrings = ['FractalDim','Fractal Dimension']
        thelabels = ['1 wall','2 walls','3 walls','4 walls']        
        
    if FlagFit : #PROIBINDO CONTINUAR IMPLEMENTAÇÃO NO CASO NÃO PREPARADO
        if kind is not 'A' : 
            print("Só é permitido FlagFit=True para kind='A'")
            return "ERROR"

    Tini = time.time()
    
    pf_range = np.linspace(0,1,prob_points)
    frequency = np.array([average_main(L,pf,samples,kind) for pf in pf_range])

    Tend = time.time()-Tini

    data_header = "Grid Size: "+str(L)+"\n"
    data_header += "Method: "+TheMethod[method]+"\n"
    data_header += "Number of samples: "+str(samples)+"\n"
    data_header += "Run time: "+str("{:.2f}".format(Tend))+" s\n"
    data_header += "Free Path Probability (P_F), "+thestrings[1]+"\n"
    data_points = np.column_stack((pf_range, frequency))
    np.savetxt(str(method)+thestrings[0]+str(L)+'.dat', data_points,header=data_header)   

    if(FlagPlot):
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        if frequency.ndim == 1 : 
            ax.plot(pf_range,frequency,marker='o',ls=':')
        elif frequency.ndim == 2 :
            for i in range(frequency.shape[1]):
                ax.plot(pf_range,frequency[:,i],label=thelabels[i])
            ax.legend(loc='upper left',fontsize='small')
        ax.set_xlabel('Free cell probability')
        ax.set_ylabel(thestrings[1]+'\n('+str(samples)+' samples)')
        ax.set_title('Percolation \n L = '+str(L)+', time = '+str("{:.2f}".format(Tend))+'s')
        fig.savefig(str(method)+thestrings[0]+str(L)+'.png')    
         
    if(FlagFit): #Funcional somente para o caso 'A'
        popt,pcov = curve_fit(sigmoid,pf_range,frequency,p0=[10,.6])
        xlin = np.linspace(0,1,200)
        fig = plt.figure()
        axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        axes1.plot(pf_range,frequency,marker='o',ls=':',label='Simulated')
        axes1.plot(xlin,sigmoid(xlin,*(popt)),'r-',label='Sigmoid\n $\sigma$[%5.2f (x-%5.5f)]'% tuple(popt))
        axes1.set_xlabel('Free cell probability')
        axes1.set_ylabel(thestrings[1]+'\n('+str(samples)+' samples)')
        axes1.set_title('Percolation\n L = '+str(L)+', time = '+str("{:.2f}".format(Tend))+'s')
        axes1.legend(loc='upper left',fontsize='small')
        fig.savefig(str(method)+'Fit'+thestrings[0]+str(L)+'.png')
    
    return Tend
   
### -------------- Function : necessário para o fit ---------------
def sigmoid(x,a,b):
    return 1 / (1+ np.exp(-a*(x-b)))


### -------------- Function : Proposta para medir dimensão fractal ---------------
import scipy.misc
def fractal_dimension(Z, threshold=0.9):
    #From https://gist.github.com/viveksck/1110dfca01e4ec2c608515f0d5a5b1d1

    # Only for 2d image
    assert(len(Z.shape) == 2)

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(Z, k):
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
                               np.arange(0, Z.shape[1], k), axis=1)

        # We count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])


    # Transform Z into a binary array
    Z = (Z < threshold)

    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))

    # Extract the exponent
    n = int(np.log(n)/np.log(2))

    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)

    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(Z, size))

    # Fit the successive log(sizes) with log (counts)
    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    return -coeffs[0]

###############################################################################

#Tentativa de simplificar codigo...

#
#SAVE_DICT = {
#        "A" : { "strings" : ['Percolation','Frequency'] },
#        "B" : { "strings" : ['Density','Mean Density'] },
#        "C" : { "strings" : ['PercWhich','Frequency'],
#               "labels" : ['up','down','left','right'] },
#        "D" : { "strings" : ['PercMany','Frequency'],
#               "labels" : ['1 wall','2 walls','3 walls','4 walls'] },
#        "E" : { "strings" : ['FractalDim','Fractal Dimension']}
#        }

#k = "C"
#if ("labels" in SAVE_DICT[k]) :
#    for value in SAVE_DICT[k]["labels"]:
#        print(value)

#IGNORAR TUDO ABAIXO...
#def SAVE_DATA():
#    data_header = "Grid Size: "+str(L)+"\n"
#    data_header += "Method: "+TheMethod[method]+"\n"
#    data_header += "Number of samples: "+str(samples)+"\n"
#    data_header += "Run time: "+str("{:.2f}".format(Tend))+" s\n"
#    data_header += "Free Path Probability (P_F), "+thestrings[1]+"\n"
#    data_points = np.column_stack((pf_range, frequency))
#    np.savetxt(str(method)+thestrings[0]+str(L)+'.dat', data_points,header=data_header)
#    return print("Data saved into file: ")
#
#def SAVE_FIT():
#    popt,pcov = curve_fit(sigmoid,pf_range,frequency,p0=[10,.6])
#    xlin = np.linspace(0,1,200)
#    fig = plt.figure()
#    axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#    axes1.plot(pf_range,frequency,marker='o',ls=':',label='Simulated')
#    axes1.plot(xlin,sigmoid(xlin,*(popt)),'r-',label='Sigmoid\n $\sigma$[%5.2f (x-%5.5f)]'% tuple(popt))
#    axes1.set_xlabel('Free cell probability')
#    axes1.set_ylabel(thestrings[1]+'\n('+str(samples)+' samples)')
#    axes1.set_title('Percolation\n L = '+str(L)+', time = '+str("{:.2f}".format(Tend)))
#    axes1.legend(loc='upper left',fontsize='small')
#    fig.savefig(str(method)+'Fit'+thestrings[0]+str(L)+'.png')
#    return
#
#def SAVE_FIG():
#    fig = plt.figure()
#    axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#    if ((kind=='A') or (kind=='B')):
#        axes1.plot(pf_range,frequency,marker='o',ls=':')
#    if ((kind=='C') or (kind=='D')):
#        axes1.plot(pf_range,frequency[:,0],label=thelabels[0])
#        axes1.plot(pf_range,frequency[:,1],label=thelabels[1])
#        axes1.plot(pf_range,frequency[:,2],label=thelabels[2])
#        axes1.plot(pf_range,frequency[:,3],label=thelabels[3])
#        axes1.legend(loc='upper left',fontsize='small')
#    axes1.set_xlabel('Free cell probability')
#    axes1.set_ylabel(thestrings[1]+'\n('+str(samples)+' samples)')
#    axes1.set_title('Percolation \n L = '+str(L)+', time = '+str("{:.2f}".format(Tend)))
#    fig.savefig(str(method)+thestrings[0]+str(L)+'.png')   
#    return
#
#def runover_pf():
#    data_output = []
#    pf_range = np.linspace(0,1,prob_points)
#    for pf in pf_range:
#        
#        data_output.append(average_main(L,pf,samples,kind))
#    return
#
#def average_pf(L,pf,samples):
#    average = 0
#    for it in range(samples):
#        x = lattice(L,pf)
#        x.percolate2()
#        average += x.edgeB
#    average /= samples
#    return average
    

##codigo desejado
#def plain(L,pf_points,samples):
#    data_output = data(L, pf in np.linspace(0,1,pf_points),samples)    
#    
#    save_data(data_output)
#    save_fig(data_output)
#    
#    return
#
#def data(L,pf,samples):
#    data=[]
#    for i in range(samples):
#        x.lattice(L,pf)
#        x.percolate()
#        data += [x.edgeB, x.density(),..]
#    data /= samples
#    return data
