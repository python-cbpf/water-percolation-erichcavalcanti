def try_time():
    L = 10
    for i in range(10000):
        N = L*L
        X = lattice(L)
        StatsNoccupied = np.zeros((N-1,3))        
        for n in range(N-1): #modifying the number of occupied cells
            X.addcell()
            StatsNoccupied[n] = np.array(X.stats())
    return

########################################################################
import time
def time_test():
    L_range = np.array([50,100,150,200,250])
    T_range = np.zeros_like(L_range,dtype=float)
    for i in range(len(L_range)):
        T_range[i] = time.time()
        one_run(L_range[i])
        T_range[i] = time.time()-T_range[i]
    return L_range,T_range

########################################################################
def test(L,many,CUT):
    """ So plota"""
    Ns,sNs,s2Ns = many_runs(L,many,CUT).T
    n = list(range(1,L*L))

    ax1 = plt.subplot(311)
    plt.plot(n,Ns,'r')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax2 = plt.subplot(312, sharex=ax1)
    plt.plot(n,sNs,'r')
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax3 = plt.subplot(313, sharex=ax1)
    plt.plot(n,s2Ns,'r')
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.show()
    return

def many_runs(L,how_many,CUT):
    """Tira a media para um certo numero de samples"""
    Stats = np.zeros((L*L-1,3))
    for runs in range(how_many):
        Stats += one_run(L,CUT)
    return Stats/how_many

def one_run(L,CUT):
    """Testa uma rede de tamanho L, todas as ocupacoes"""
    N = L*L
    X = lattice(L)
    StatsNoccupied = np.zeros((N-1,3))
    
    for n in range(N-1): #modifying the number of occupied cells
        X.addcell()
        StatsNoccupied[n] = np.array(X.stats_cut()[CUT])
    
    return StatsNoccupied

########################################################################
import numpy as np
class lattice:
    def __init__(self,L):
        """ 
        L, grid size
        Initialize empty grid and empty cluster dictionary
        """
        self.L = L  
        self.GRID = np.zeros((self.L,self.L),dtype=bool)
        self.CLUSTERS = {}
        self.NEXT = 1 #Takes control of CLUSTER label
        self.listindex() #Guarantee that it is called
        return
    def listindex(self):
        """ Get a list of cells in the order to be added """
        self.ALLINDEX = list(zip(*np.where(self.GRID==0))) ##otimização. ZIP junta o obtido com where
        np.random.shuffle(self.ALLINDEX)
        return
    def addcell(self):
        if len(self.ALLINDEX)==1: return

        INDEX = self.ALLINDEX.pop(0) #Get one index from the list
        self.GRID[INDEX] = True #Activate cell in the grid
        
        #confere a que cluster pertencem os vizinhos
        #Necessario tomar cuidado para não contar o mesmo cluster varias vezes
        WHICH_CLUSTER = []
        ALL_NGH = self.neighbor_D(INDEX)
        for NGH in ALL_NGH:
            AUX = self.belongsto(NGH) 
            if len(AUX)==1:
                if (AUX[0] not in WHICH_CLUSTER): WHICH_CLUSTER += AUX
        
        if(len(WHICH_CLUSTER)==0): #nova celula não tem clusters vizinhos
            self.CLUSTERS[self.NEXT] = [INDEX] #list of tuples
            self.NEXT+=1 #new cluster created, therefore INDEX is incremented
        elif(len(WHICH_CLUSTER)==1): #apenas um vizinho
            self.CLUSTERS[WHICH_CLUSTER[0]] += [INDEX]
        else : #mais de um vizinho
            self.CLUSTERS[WHICH_CLUSTER[0]] += [INDEX] #arbitrarily add to first cluster
            for EACH in WHICH_CLUSTER[1:]:
                self.CLUSTERS[WHICH_CLUSTER[0]] += self.CLUSTERS.pop(EACH)
        return
    def belongsto(self,INDEX):
        """ verifica a que cluster pertence um elemento. """
        if self.GRID[INDEX] : #confere se celula esta ativa
            for LABEL in self.CLUSTERS:
                if INDEX in self.CLUSTERS[LABEL]:
                    return [LABEL]
        return []
    def neighbor_D(self,INDEX):
        """ Boundary empty"""        
        NGH = []        
        if INDEX[0]!=0 : NGH += [(INDEX[0]-1,INDEX[1])] 
        if INDEX[0]!=self.L-1 : NGH += [(INDEX[0]+1,INDEX[1])]             
        if INDEX[1]!=0 : NGH += [(INDEX[0],INDEX[1]-1)] 
        if INDEX[1]!=self.L-1 : NGH += [(INDEX[0],INDEX[1]+1)] 
        return NGH
    def neighbor_P(self,INDEX):
        """ Periodic boundary condition"""
        if INDEX[0]==0 : NGH = [(self.L-1,INDEX[1])]
        else : NGH = [(INDEX[0]-1,INDEX[1])] 
    
        if INDEX[0]==self.L-1 : NGH += [(0,INDEX[1])]
        else : NGH += [(INDEX[0]+1,INDEX[1])] 
            
        if INDEX[1]==0 : NGH += [(INDEX[0],self.L-1)]
        else : NGH += [(INDEX[0],INDEX[1]-1)] 

        if INDEX[1]==self.L-1 : NGH += [(INDEX[0],0)]
        else : NGH += [(INDEX[0],INDEX[1]+1)] 

        return NGH
    def check_InfClust(self):
        """ Check if some cluster touch itself through edges

            - Feio do modo como esta
        """
        #Modo nao inteligente de construir pares e conferir se existe cluster infinito
        pairs = [[(0,i),(self.L-1,i)] for i in range(self.L)] + [[(i,0),(i,self.L-1)] for i in range(self.L)]
        
        InfClust = []
        if len(pairs)!=0:
            for eachclust in self.CLUSTERS:
                if (len(self.CLUSTERS[eachclust])>(self.L-1)): #essencial para nao olhar clusters pequenos
                    for eachpair in pairs:
                        flag = (eachpair[0] in self.CLUSTERS[eachclust]) & (eachpair[1] in self.CLUSTERS[eachclust])
                        if flag : InfClust += [eachclust]; break
        return InfClust
    def stats_cut(self):
        newdict = {}
        
        for label in self.CLUSTERS:
            newdict[label] = len(self.CLUSTERS[label])
        newdict_cut = newdict.copy()
        for label in self.check_InfClust():
            newdict_cut.pop(label)
            
        LISTCLST, LSITCLST_CUT = np.array(list(newdict.values())), np.array(list(newdict_cut.values()))
        STATS = [len(LISTCLST), LISTCLST.sum(), (LISTCLST*LISTCLST).sum()]
        STATS_CUT = [len(LSITCLST_CUT), LSITCLST_CUT.sum(), (LSITCLST_CUT*LSITCLST_CUT).sum()]
        return STATS,STATS_CUT
    def stats(self): ##Parece ser mais veloz que o stats
        newdict = {}
        
        for label in self.CLUSTERS:
            newdict[label] = len(self.CLUSTERS[label])
            
        LISTCLST = np.array(list(newdict.values()))
        STATS = [len(LISTCLST), LISTCLST.sum(), (LISTCLST*LISTCLST).sum()]
        return STATS
from matplotlib import pyplot as plt
def plot(toprint):
    plt.figure()
    plt.imshow(toprint,vmin=-1,vmax=1,cmap='coolwarm')
    plt.axis('off')
    return

########################################################################
