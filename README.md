# All

There are 3 codes at the moment

"Perc_Erich" 
    The first one (and is better "documented"). Works with the whole matrix, get the clusters from it.

"FindUnion"
    Add one cell at a time. Obtain clusters.
    
"RootPointer"
    Same as above, but with the root for each cluster


# water-percolation-erichcavalcanti
water-percolation-erichcavalcanti created by GitHub Classroom

@author: erichcavalcanti@gmail.com

Temos alguns diferentes tipos de objetivos (kind)...

    A - frequencia com que atinge no mínimo uma vez uma das bordas da rede
    B - densidade de ocupação quando ocorre a condição acima
    C - frequencia com que atinge um lado especifico da rede [up,down,left,right]
    D - frequencia com que atinge um número especifico de lados da rede [1,2,3,4]
    E - "dimensão fractal"

Processo: runover_main > average_main > chama classe > nextstep, check

    Função nextstep responsável pela evolução da rede
    Função check confere se evolução deve parar
    A classe gera a rede inicial e evolui ela. Gera a figura da rede final.
    Função average_main repete main um certo número de samples para tirar
        média de alguma quantidade em Determinada frequencia fixa (ver objetivos)
    Função runover_main pega average_main e varia a frequencia.
    
    ----------------------------------------------------------------
após compilar funções, rodar:

    runover_main(L,samples,prob_points,kind,FlagFit,FlagPlot)
exemplos:

    runover_main(21,1000,10,'A',True,True)
    runover_main(21,1000,10,'D',False,True)

OU de modo explicito:

    x.lattice(L,p) #gera rede original
    x.percolate() #roda o código
    plot(x.toprint()) #mostra figura
    x.edgeB #diz se atingiu ao menos uma borda
    x.edge #diz quandas bordas atingiu
    x.density() #diz densidade
    fractal_dimension(x.toprint(),.9) #mede (?) dimensão fractal (?) da figura
    
---------------------------------------------------------------------
To-Do-List
[ ] Medir Cluster-size

Mais implementações possíveis:
[ ] Redes com conctividades diferentes
[ ] Finite-size scaling. Observar cruzamento das curvas
[ ] Adicionar incertezas
[ ] Desenhar rede fractal
[ ] Implementar pandas para lidar com todos os possíveis outputs de lattice?


https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.51.2347
https://pdfs.semanticscholar.org/707e/2fe48cd821117c5e142d070a164b0b8f6156.pdf
https://arxiv.org/pdf/1504.02898.pdf

