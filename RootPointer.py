########################################################################
import numpy as np
class lattice:
    def __init__(self,L):
        """ 
        L, grid size
        Initialize: empty grid, grid of pointers, list of index
        """
        self.L = L  
        self.GRID = np.zeros((self.L,self.L),dtype=bool)
        self.DIR_GRID = np.zeros((self.L,self.L),dtype=np.complex_)
        self.PERC_GRID = self.GRID.copy()
        self.POS = np.mgrid[0:self.L,0:self.L][1]-np.mgrid[0:self.L,0:self.L][0]*1j
        self.listindex()
        return
    def listindex(self):
        """ Get a list of cells in the order to be added """
        self.ALLINDEX = list(zip(*np.where(self.GRID==0)))
        np.random.shuffle(self.ALLINDEX)
        return
    def addcell(self):
        """ Remember that Y axis has an inverted direction 
        """
        if len(self.ALLINDEX)==1: return
        
        #Get one index from listindex and activate cell in the grid
        INDEX = self.ALLINDEX.pop(0)
        self.GRID[INDEX] = True #root of itself
    
        # Neighbors of the tested cell
        ALL_NGH_ROOT, DUPLICATED_ROOT = self.checkneighbor(INDEX)
        
        # New pointers
        LEN = len(ALL_NGH_ROOT)
        if LEN >0:
            #Points to the first option
            self.attachto(ALL_NGH_ROOT[0],INDEX)
            #Now make all other roots point to the new root
            for NGH_ROOT in ALL_NGH_ROOT[1:]:
                self.attachto(ALL_NGH_ROOT[0],NGH_ROOT)
        #refresh grid IF we have a union between two different clusters
        if LEN>1 : self.refresh()
        
        # Percolation Matrix (tells if cell belongs to a percolated cluster)  
        NBH_HAS_PERC = [self.PERC_GRID[self.BC_P(NGH_ROOT)] for NGH_ROOT in ALL_NGH_ROOT]
        if(len(NBH_HAS_PERC)>0):
            # if all neighbors are percolated clusters
            if all(NBH_HAS_PERC):
                self.PERC_GRID[INDEX] = True
            # if at least one neighbor is a percolated cluster
            elif LEN>1:
                if any(NBH_HAS_PERC):
                    self.PERC_GRID[INDEX] = True
                    for NGH_ROOT in ALL_NGH_ROOT : self.PERC_GRID[self.BC_P(NGH_ROOT)] = True
                    self.refresh_perc()
            # "union" between the same cluster. Check percolation.
            elif len(DUPLICATED_ROOT)>0 : 
                if self.ACT_NBH_DIFF(INDEX): #Check percolation
                    self.PERC_GRID[INDEX] = True
                    for NGH_ROOT in ALL_NGH_ROOT : self.PERC_GRID[self.BC_P(NGH_ROOT)] = True
                    self.refresh_perc()
        return
    def ACT_NBH_DIFF(self,INDEX):
        # confere se diferença de distancia para o qualquer um dos primeiros vizinhos é maior do que 1
        # essa diferença demonstra diferença de caminho até a raiz e, portanto, percolação
        AUX = [(INDEX[0]-1,INDEX[1]),(INDEX[0]+1,INDEX[1]),(INDEX[0],INDEX[1]-1),(INDEX[0],INDEX[1]+1)]
        DIST = []
        for EACH in AUX:
            if self.GRID[self.BC_P(EACH)]:
                DIST += [abs(self.DIR_GRID[self.BC_P(EACH)]-self.DIR_GRID[INDEX])]
        
        return any(np.array(DIST)>1.)
    def checkneighbor(self,INDEX):
        #List of ROOTS of the active neighbors
        ALL_NGH_ROOT = self.neighbor_root(INDEX)#).sort(reverse=True) #sort para codigo abaixo
        DUPLICATED_ROOT = []

        if len(ALL_NGH_ROOT)>0:
            #Remove duplicates      
            AUX_LIST = [(int(ELEM[0]%self.L),int(ELEM[1]%self.L)) for ELEM in ALL_NGH_ROOT]
    
            for ELEM in AUX_LIST:
                while AUX_LIST.count(ELEM)>1 :
                   ALL_NGH_ROOT.pop( AUX_LIST.index(ELEM) )
                   AUX_LIST.remove(ELEM)
                   if ELEM not in DUPLICATED_ROOT : DUPLICATED_ROOT += [ELEM]
                   
        return ALL_NGH_ROOT,DUPLICATED_ROOT
    def attachto(self,ROOT,INDEX):
        #Points to the first option
        self.DIR_GRID[self.BC_P(INDEX)] = (ROOT[1]-INDEX[1]) - (ROOT[0]-INDEX[0])*1j
        return
    def BC_P(self,INDEX): #condicao periodica
        return (int(INDEX[0]%self.L),int(INDEX[1]%self.L))  #tuple( np.array(INDEX) %L)
    def NEXTCELL(self,INDEX):
        return (int(INDEX[0]-self.DIR_GRID[self.BC_P(INDEX)].imag),int(INDEX[1]+self.DIR_GRID[self.BC_P(INDEX)].real))
    def refresh_perc(self):
        #confere qual a matriz de "roots"
        auxind = (self.DIR_GRID + self.POS)*self.GRID
        #obriga a ter mesmo valor que seu ROOT (se ativo)
        self.PERC_GRID = self.GRID*self.PERC_GRID[np.array(-auxind.imag,dtype='i8')%self.L,np.array(auxind.real,dtype='i8')%self.L]
        return
    def refresh(self):
        #confere qual a matriz de "roots" (falsos roots na verdade)
        auxind = (self.DIR_GRID + self.POS)*self.GRID
        #obtem o pointer dos roots e adiciona no pointer presente
        self.DIR_GRID += self.GRID*(self.DIR_GRID!=0)*self.DIR_GRID[np.array(-auxind.imag,dtype='i8')%self.L,np.array(auxind.real,dtype='i8')%self.L]
        return
    def neighbor(self,INDEX): #only active ones
        NGH = [] 

        AUX_NGH = (INDEX[0]-1,INDEX[1]) #Indice do vizinho
        if self.GRID[self.BC_P(AUX_NGH)] : NGH += [AUX_NGH]

        AUX_NGH = (INDEX[0]+1,INDEX[1])             
        if self.GRID[self.BC_P(AUX_NGH)] : NGH += [AUX_NGH]

        AUX_NGH = (INDEX[0],INDEX[1]-1) 
        if self.GRID[self.BC_P(AUX_NGH)] : NGH += [AUX_NGH]

        AUX_NGH = (INDEX[0],INDEX[1]+1) 
        if self.GRID[self.BC_P(AUX_NGH)] : NGH += [AUX_NGH]

        return NGH
    def neighbor_root(self,INDEX): #root of activated neighbors
        NB_R = [] 

        AUX = (INDEX[0]-1,INDEX[1]) #Indice do vizinho
        if self.GRID[self.BC_P(AUX)] : NB_R += [self.NEXTCELL(AUX)]

        AUX = (INDEX[0]+1,INDEX[1])             
        if self.GRID[self.BC_P(AUX)] : NB_R += [self.NEXTCELL(AUX)]

        AUX = (INDEX[0],INDEX[1]-1) 
        if self.GRID[self.BC_P(AUX)] : NB_R += [self.NEXTCELL(AUX)]

        AUX = (INDEX[0],INDEX[1]+1) 
        if self.GRID[self.BC_P(AUX)] : NB_R += [self.NEXTCELL(AUX)]

        return NB_R
    def cluster_find(self): #lista de indices de todos os clusters
        return list(zip(*np.where(self.GRID*(self.DIR_GRID==0))))
    def cluster_info(self):
        #lista de roots para cada cluster
        CLUSTER_LIST = self.cluster_find()
        #grid informando root, cuidado pois celula vazia usa (0,0)
        AUX_GRID = (self.DIR_GRID + self.POS)*self.GRID
        CLUSTER_GRID = (AUX_GRID.real%self.L) + (AUX_GRID.imag%self.L) *1j
        
        #Obtem tamanho de cada cluster. 
        CLUSTER_SIZE = []
        CLUSTER_PERC = []
        for EACH_CLUSTER in CLUSTER_LIST:
            #conversao de tuple para complexo
            value = (EACH_CLUSTER[1]%self.L) + ((-EACH_CLUSTER[0])%self.L)*1j
            if value !=0 : 
                CLUSTER_SIZE += [(CLUSTER_GRID==value).sum()]
            else : 
                CLUSTER_SIZE += [(CLUSTER_GRID==value).sum()-(self.GRID==value).sum()]
            # aproveitando para guardar informação se percolou ou não
            CLUSTER_PERC += [self.PERC_GRID[self.BC_P(EACH_CLUSTER)]]
        return [CLUSTER_LIST, CLUSTER_SIZE, CLUSTER_PERC]
    def stats(self):
        CLUSTER_LIST, CLUSTER_SIZE, CLUSTER_PERC = self.cluster_info()
        SIZE_ARRAY = np.array(CLUSTER_SIZE)
        return [len(SIZE_ARRAY), SIZE_ARRAY.sum(), (SIZE_ARRAY*SIZE_ARRAY).sum()]
    def stats_cut(self):
        CLUSTER_LIST, CLUSTER_SIZE, CLUSTER_PERC = self.cluster_info()
        
        #Aplicando o CUT
        i=0
        for elem in list(*np.where(np.array(CLUSTER_PERC))):
            CLUSTER_SIZE.pop(elem-i)
            i+=1
        
        SIZE_ARRAY = np.array(CLUSTER_SIZE)
        return [len(SIZE_ARRAY), SIZE_ARRAY.sum(), (SIZE_ARRAY*SIZE_ARRAY).sum()]
  
########################################################################
from matplotlib import pyplot as plt
def plot(toprint):
    plt.figure()
    plt.imshow(toprint,vmin=-1,vmax=1,cmap='coolwarm')
    plt.axis('off')
    return
def Vplot(grid,vecgrid):
    plt.figure()
    plt.imshow(grid,vmin=0,vmax=2,cmap='coolwarm')
    plt.quiver(np.real(vecgrid),np.imag(vecgrid),color='green',units='xy')
    plt.axis('off')
    return

########################################################################
def DOIT_percolate(L):
    X = lattice(L)
    for i in range(L*L):
        X.addcell()
        if X.PERC_GRID.any() : 
            print(i)
            Vplot(X.GRID*1+((X.PERC_GRID==1)*X.GRID)*(1),X.DIR_GRID)
            Vplot(X.GRID*1+((X.DIR_GRID==0)*X.GRID)*(1),X.DIR_GRID)
            break
    return
########################################################################

def this_onerun(L):
    N = L*L
    X = lattice(L)
    StatsNoccupied = np.zeros((N-1,3))
    
    for n in range(N-1): #modifying the number of occupied cells
        X.addcell()
        StatsNoccupied[n] = np.array(X.stats)
    
    return StatsNoccupied
