from sage.all import Graph, ZZ, Combinations, prod, factorial
import itertools
import time
from .utils import float2time

class LabeledStableGraph:
    '''
    Represents Labeled Stable Graph -- the canonical representative of the labeled stable graph given by edges, loops and kappa.
    
    INPUT:

        - ``edges``  -- tuple of triples, a triple (v1,v2,m) means that the vertices v1 and v2 are connected by m edges
        - ``loops``  -- list, an integer loops[i] is the number of loops associated to the vertex i
        - ``kappa``  -- list of lists, a partition of stratum into sublists, where a list kappa[i] is a list of orders of zeroes           associated to the vertex i
    '''
    def __init__(self,edges,loops,kappa):
        if not edges:
            graph = Graph(weighted=True, loops=False, multiedges=False)
            graph.add_vertex()
        else:
            graph = Graph(list(edges), loops=False, multiedges=False, weighted=True)        
        self.edges, self.loops, self.kappa, self.graph = self.canonical(edges,loops,kappa,graph)
        self.genera = [(sum(self.kappa[v])-2*self.vertex_deg(v)-4*self.loops[v]+4)/ZZ(4) for v in self.graph.vertices()]
    
    def k_to_p(self,edges,loops,kappa,graph):
        '''
        Return a canonical partition of vertices of the graph into lists grouped by the same number of loops and orders of               zeroes. The order is lexicographical with respect to kappa -> loops -> edges.
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
    
    def canonical(self,edges,loops,kappa,graph):
        '''
        Return a 4-tuple (edges,loops,kappa,graph), corresponding to the canonical representative of the class of isomorphism of 
        the given labeled stable graph, where only vertices with the same number of loops and zero orders are allowed to permute         and only edges of the same weight are allowed to permute.
        '''
        can_gr, relab = graph.canonical_label(partition=self.k_to_p(edges,loops,kappa,graph), certificate=True, edge_labels=True)
        can_loops = list(loops)
        can_kappa = list(kappa)
        for k,v in relab.items():
            can_loops[v] = loops[k] 
            can_kappa[v] = kappa[k]
        return tuple(can_gr.edges()), tuple(can_loops), tuple(can_kappa), can_gr.copy(immutable=True)
   
    def __repr__(self):
        return f"Labeled Stable Graph ({self.edges},{self.loops},{self.kappa})"
        
    def vertex_deg(self, v):
        '''
        Return the sum of the weights at the vertex v of the graph given by the list of edges.
        '''
        deg = 0
        for e in self.edges:
            if e[0]==v or e[1]==v:
                deg += e[2]
        return deg
        
    def genera(self):
        '''
        Return a list of genera of vertices of the labeled stable graph.
        '''
        return [(sum(self.kappa[v])-2*self.vertex_deg(v)-4*self.loops[v]+4)/ZZ(4) for v in self.graph.vertices()]
    
    def loop_degenerations(self):
        '''
        Return the set of all Labeled Stable Graphs, obtained by the degenerations adding a loop to this Labeled Stable Graph.
        '''
        new_graphs = set()
        for v in self.graph.vertices():
            assert self.genera[v] in ZZ, f"The genus of the {repr(self)} at vertex {v} is not an integer."
            if self.genera[v]>0:
                new_loops = list(self.loops)
                new_loops[v] = self.loops[v]+1
                new_graphs.add(LabeledStableGraph(self.edges,new_loops,self.kappa))
        return new_graphs

    def edge_degenerations(self):
        '''
        Return the set of all Labeled Stable Graphs, obtained by *special* degenerations adding an edge to this Labeled Stable           Graph. A *special* edge degeneration only adds an edge to the vertex with largest profile, when it's unique. Profile is           the sum of weight and length of the associated zeros partition.
        '''
        new_graphs = set()
        profile = [sum(self.kappa[v])+len(self.kappa[v]) for v in self.graph.vertices()]
        m = max(profile)
        if profile.count(m) == 1: # add an edge only to the unique vertex with 'maximal' zeroes partition
            v = profile.index(m)
            v_edges = self.graph.edges_incident(v)
            v_edge_weights = [i[2] for i in v_edges]
            v_comb_weights = [[[weight-i,i] for i in range(weight+1)] for weight in v_edge_weights]
            if not v_comb_weights: v_Split_Weights = [()]
            v_Split_Weights = list(itertools.product(*v_comb_weights))
            v_Split_Loops = [[self.loops[v]-i-j,j,i] for i in range(self.loops[v]+1) for j in range(self.loops[v]+1-i)]
            if not v_Split_Loops: v_Split_Loops = [[0,0,0]]
            v_Subset_Zeroes = [i for i in Combinations(self.kappa[v]) if sum(i)%2==0][1:-1]
            for split_weights,split_loops,sub_zeroes in itertools.product(v_Split_Weights,v_Split_Loops,v_Subset_Zeroes):
                new_genus_0 = (sum(sub_zeroes)-2*sum(i[0] for i in split_weights)-2*(split_loops[2]+1)-4*split_loops[0]+4)/ZZ(4)
                new_genus_1 = (sum(self.kappa[v])-sum(sub_zeroes)-2*sum(i[1] for i in split_weights)-2*(split_loops[2]+1)-4*split_loops[1]+4)/ZZ(4)
                if self.genera[v]>=new_genus_0>=0 and self.genera[v]>=new_genus_1>=0 and new_genus_0 in ZZ and new_genus_1 in ZZ:
                    # creating copies
                    new_graph = self.graph.copy(immutable=False)
                    new_loops = list(self.loops)
                    new_kappa = list(self.kappa)
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
                    new_graphs.add(LabeledStableGraph(new_graph.edges(),new_loops,new_kappa))
        return new_graphs


    def one_step_degenerations(self):
        '''
        Return the set of all Labeled Stable Graphs, obtained by *special* one step degenerations of this Labeled Stable Graph. A         *special* degeneration is adding a loop to a single vertex graph, or adding an edge to any graph.
        '''
        degenerations = set()
        if len(self.loops) == 1:    # make loop degeneration only if the graph has a single vertex
            degenerations.update(self.loop_degenerations())
        degenerations.update(self.edge_degenerations())
        return degenerations
    
    def Aut(self):
        '''
        Return the order of the group of automorphisms of the labeled stable graph.
        '''
        graph_aut = self.graph.automorphism_group(partition=self.k_to_p(),edge_labels=True,order=True,return_group=False)
        loops_aut = prod(2**i*factorial(i) for i in self.loops)
        weights_aut = prod(factorial(e[2]) for e in self.edges)
        return graph_aut*loops_aut*weights_aut

def stable_lab_graphs(stratum, by_codim=False, one_vertex=False, verbose=False):
    '''
    Return a list of all labeled stable graphs given by the stratum.
    '''
    assert sum(stratum)%4 == 0, f"The sum of orders of zeroes of the stratum has to be a multiple of 4."
    kappa = [sorted(stratum)]
    codim = 0
    g = (sum(stratum)+4)/4
    
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
    degenerations = [ {LabeledStableGraph([],[0],kappa)} ]
    while degenerations[codim]:
        tic = time.time()
        codim += 1
        # adding empty list for degenerations of the next codimension 
        degenerations.append(set())
        # looping though graphs of previous codimension
        for stg in degenerations[codim-1]:
            # taking a union with all possible one-step degenerations of each graph
            degenerations[codim].update(stg.one_step_degenerations())
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
