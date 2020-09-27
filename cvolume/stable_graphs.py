from sage.all import Graph, ZZ, Combinations, prod, factorial
import itertools
import time
from .utils import float2time

def init_graph(edges):
    '''
    Initialize a graph from the list of edges. Return both edges and a graph.
    '''
    if not edges:
        graph = Graph(weighted=True, loops=False, multiedges=False)
        graph.add_vertex()
    else:
        graph = Graph(list(edges), loops=False, multiedges=False, weighted=True)
    return list(graph.edges()), graph

def vertex_deg(edges,v):
    '''
    Return the sum of the weights at the vertex v of the graph given by the list of edges.
    '''
    deg = 0
    for e in edges:
        if e[0]==v or e[1]==v:
            deg += e[2]
    return deg

def genera(edges,loops,kappa,graph):
    '''
    Return a list of genera of vertices of the labeled stable graph (edges,loops,kappa,graph).
    '''
    return [(sum(kappa[v])-2*vertex_deg(edges,v)-4*loops[v]+4)/ZZ(4) for v in graph.vertices()]

def add_loop(edges,loops,kappa,graph):
    '''
    Return a list of all labeled stable graphs, obtained by adding a loop to the labeled stable graph (edges, loops, kappa, graph).
    '''
    new_graphs = []
    #genus = [(sum(kappa[v])-2*vertex_deg(edges,v)-4*loops[v]+4)/4 for v in graph.vertices()]
    genus = genera(edges,loops,kappa,graph)
    for v in graph.vertices():
        assert genus[v] in ZZ, f"The genus of graph {(graph.edges(),loops,kappa)} at vertex {v} is not an integer."
        if genus[v]>0:
            new_loops = list(loops)
            new_loops[v] = loops[v]+1
            new_graphs.append((edges,tuple(new_loops),kappa,graph))
    return new_graphs

def add_edge(edges,loops,kappa,graph):
    '''
    Return a list of all labeled stable graphs, obtained by adding an edge to the labeled stable graph (edges, loops, kappa, graph).
    '''
    new_graphs = []
    genus = genera(edges,loops,kappa,graph)
    profile = [sum(kappa[v])+len(kappa[v]) for v in graph.vertices()]
    m = max(profile)
    if profile.count(m) == 1: 
        v = profile.index(m)
        v_edges = graph.edges_incident(v)
        v_edge_weights = [i[2] for i in v_edges]
        v_comb_weights = [[[weight-i,i] for i in range(weight+1)] for weight in v_edge_weights]
        if not v_comb_weights: v_Split_Weights = [()]
        v_Split_Weights = list(itertools.product(*v_comb_weights))
        v_Split_Loops = [[loops[v]-i-j,j,i] for i in range(loops[v]+1) for j in range(loops[v]+1-i)]
        if not v_Split_Loops: v_Split_Loops = [[0,0,0]]
        v_Subset_Zeroes = [i for i in Combinations(kappa[v]) if sum(i)%2==0][1:-1]
        for split_weights,split_loops,sub_zeroes in itertools.product(v_Split_Weights,v_Split_Loops,v_Subset_Zeroes):
            new_genus_0 = (sum(sub_zeroes)-2*sum(i[0] for i in split_weights)-2*(split_loops[2]+1)-4*split_loops[0]+4)/ZZ(4)
            new_genus_1 = (sum(kappa[v])-sum(sub_zeroes)-2*sum(i[1] for i in split_weights)-2*(split_loops[2]+1)-4*split_loops[1]+4)/ZZ(4)
            if genus[v]>=new_genus_0>=0 and genus[v]>=new_genus_1>=0 and new_genus_0 in ZZ and new_genus_1 in ZZ:
                # creating copies
                new_graph = graph.copy(immutable=False)
                new_loops = list(loops)
                new_kappa = list(kappa)
                new_kappa[v] = list(new_kappa[v])
                # add new vertex
                new_graph.add_vertex()
                new_v = new_graph.vertices()[-1]
                # distribute zeroes
                for i in sub_zeroes:
                    new_kappa[v].remove(i)
                new_kappa.append(tuple(sorted(new_kappa[v])))
                new_kappa[v] = tuple(sorted(sub_zeroes))              
                # distribute loops and add an extra edge
                new_loops[v] = split_loops[0]
                new_loops.append(split_loops[1]) 
                new_graph.add_edge(v,new_v,split_loops[2]+1)
                # distribute edges
                for i in range(len(v_edges)):
                    if split_weights[i][0] != 0:
                        new_graph.set_edge_label(v_edges[i][0],v_edges[i][1],split_weights[i][0])
                    else:
                        new_graph.delete_edge(v_edges[i])
                    if split_weights[i][1] != 0:
                        if v_edges[i][0] == v:
                            new_graph.add_edge(new_v,v_edges[i][1],split_weights[i][1])
                        elif v_edges[i][1] == v:
                            new_graph.add_edge(v_edges[i][0],new_v,split_weights[i][1])
                        else:
                            print(f"Edge {v_edges[i]} is not adjacent to vertex {v}")
                # add new stable graph to the list
                new_graphs.append((tuple(new_graph.edges()),tuple(new_loops),tuple(new_kappa),new_graph.copy(immutable=True)))
    return new_graphs

def k_to_p(edges,loops,kappa,graph):
    '''
    Return a canonical partition of vertices of the graph into lists grouped by the same number of loops and zeroes. 
    The order is lexicographical with respect to [kappa,loops,edges]
    '''
    edge_profile = []
    for v in graph.vertices():
        multi = [v_edge[2] for v_edge in graph.edges_incident(v)]
        edge_profile.append(sorted(multi))        
    klev = [[kappa[v],loops[v],edge_profile[v],v] for v in range(len(kappa))]
    klev = sorted(klev)
    partition = [[klev[0][3]]]
    for i in range(1,len(klev)):
        if [klev[i][0],klev[i][1],klev[i][2]] == [klev[i-1][0],klev[i-1][1],klev[i-1][2]]:
            last = partition.pop()
            last.append(klev[i][3])
            partition.append(last)
        else:
            partition.append([klev[i][3]])
    return partition

def canonical_stable(edges,loops,kappa,graph):
    '''
    Return a labeled stable graph, which is the canonical representative of the class of isomorphism of 
    the labeled stable graph (edges,loops,kappa,graph), where only vertices with the same number of loops
    and zero orders are allowed to permute and only edges of the same weight are allowed to permute.
    '''
    can_gr, relab = graph.canonical_label(partition=k_to_p(edges,loops,kappa,graph), certificate=True, edge_labels=True)
    can_loops = list(loops)
    can_kappa = list(kappa)
    for k,v in relab.items():
        can_loops[v] = loops[k] 
        can_kappa[v] = kappa[k]
    return tuple(can_gr.edges()), tuple(can_loops), tuple(can_kappa), can_gr.copy(immutable=True)

def degeneration_step(edges,loops,kappa,graph):
    '''
    Returns a list of all canonical representatives of labeled stable graphs obtained
    by all one step degenerations (adding a loop or adding an edge) of the labeled stable graph given by (edges,loops,kappa).
    '''
    degenerations = set()
    if len(loops) == 1:    # add loops only if the graph has a single vertex
        for stg in add_loop(edges,loops,kappa,graph):
            can_stg = canonical_stable(*stg)
            degenerations.add(can_stg)
    for stg in add_edge(edges,loops,kappa,graph):
        can_stg = canonical_stable(*stg)
        degenerations.add(can_stg)
#     assert not any(g[1].allows_loops() for g in degenerations)
#     assert not any(g[1].allows_multiple_edges() for g in degenerations)
    return degenerations

def stable_lab_graphs(stratum, by_codim=False, one_vertex=False, verbose=False):
    '''
    Return a list of all labeled stable graphs given by the stratum.
    '''
    assert sum(stratum)%4 == 0, f"The sum of orders of zeroes of the stratum has to be a multiple of 4."
    kappa = [tuple(sorted(stratum))]
    codim = 0
    g = (sum(stratum)+4)/4
    edges, graph = init_graph([])
    # degeneration is list whose ith element is a list of stable graphs of codimension i
    # start with a graph represented by a single vertex of genus g   
    
    ## Fix to just change the loops
#     if one_vertex:
#         degenerations = {((), (0,), tuple(kappa), graph.copy(immutable=True))}
#         for i in range(g):
#             edges,loops,kappa,graph = list(degenerations)[-1]
#             next_stg = add_loop(edges,loops,kappa,graph)[0]
#             degenerations.append(next_stg)
#         return degenerations[1:]
    degenerations = [ {((), (0,), tuple(kappa), graph.copy(immutable=True))} ]
    while degenerations[codim]:
        tic = time.time()
        codim += 1
        # adding empty list for degenerations of the next codimension 
        degenerations.append(set())
        # looping though graphs of previous codimension
        for stg in degenerations[codim-1]:
            # taking a union with all possible one-step degenerations of each graph
            degenerations[codim].update(degeneration_step(*stg))
        toc = time.time()
        if verbose: print(f"Generated {len(degenerations[-1])} codimension {codim} graphs in {float2time(toc-tic,5)}")
    if verbose:
        print(f"The total number of stable graphs for stratum {stratum} is: {sum(len(i) for i in degenerations)}.")
        toc = time.time()
        print(f"Generated all stable graphs for stratum {stratum} in: {float2time(toc-tic,5)}")  
    if by_codim:
        return degenerations[:-1] # remove the last empty list and return
    else:
        return set().union(*degenerations[1:]) # take union and return

def Aut(edges,loops,kappa):
    '''
    Return the order of the group of automorphisms of the labeled stable graph.
    '''
    edges, graph = init_graph(edges)
    graph_aut = graph.automorphism_group(partition=k_to_p(edges,loops,kappa,graph), \
                                         edge_labels=True,order=True,return_group=False)
    loops_aut = prod(2**i*factorial(i) for i in loops)
    weights_aut = prod(factorial(e[2]) for e in edges)
    return graph_aut*loops_aut*weights_aut
