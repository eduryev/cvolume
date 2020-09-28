from sage.all import factorial, Rational, prod, Combinations, zeta, pi
import itertools
import time
from .localpoly import Nlocal
from .utils import *
from .stable_graphs import *
from .series import Fs

def replmon(n):
    '''
    Return the result of application of Z-operator to the monomial b^n.
    '''
    if n==0: return 1
    else: return factorial(n)*Rational(zeta(n+1)/pi**(n+1))
    
def operator(Poly):
    '''
    Return the result of application of Z-operator to the polynomial.
    '''
    return sum(Poly.monomial_coefficient(monom)*prod(replmon(k) for k in monom.exponents()[0]) for monom in Poly.monomials())     

def graph_poly(stg):
    '''
    Return the 'Kontsevich polynomial' associated to the Labeled Stable Graph.
    
    EXAMPLES:
    
    Here we compute the polynomial associated with a labeled stable graph with one loop::
    
    sage: from cvolume import LabeledStableGraph
    sage: from cvolume.cvolume import graph_poly
    sage: stg = LabeledStableGraph([], [1], [[3, 3, -1, -1]])
    sage: graph_poly(stg)
    19/128*b1^5
    
    Here is another example for a graph with two vertices::
    
    sage: stg = LabeledStableGraph([(0,1,1)], [1, 1], [[3, -1], [3, -1]])
    sage: graph_poly(stg)
    9/16*b1*b2*b3
    
    Note that the previous example can be expressed through local polynomials associated to vertices::
    
    sage: from cvolume import Nlocal
    sage: S = PolynomialRing(QQ,['b%d' % i for i in range(1,4)])
    sage: b1,b2,b3 = S.gens()
    sage: graph_poly(stg) == 1/2*1/8*b1*b2*b3*Nlocal(0,3,[3,-1])(b1=b1,b2=b1,b3=b2)*Nlocal(0,3,[3,-1])(b1=b2,b2=b3,b3=b3)
    True
    '''
    edges,loops,kappa,graph = stg.edges,stg.loops,stg.kappa,stg.graph
    c = ZZ(1)/2**(len(graph.vertices())-1)*1/ZZ(stg.Aut())
    variables = list(S.gens())
    b = list(S.gens())
    valency = [stg.vertex_deg(v)+2*loops[v] for v in graph.vertices()]
    used_vars = []
    # dictionary that converts edges to variables
    edge_to_var = []
    for edge in edges:
        edge_vars = []
        for i in range(edge[2]):
            edge_vars.append(b[0])
            b.remove(b[0])
        edge_to_var.append((edge,edge_vars))
    edge_to_var = dict(edge_to_var)
    # list of variables assigned to loops
    loop_to_var = []
    for i in loops:
        loop_vars = []
        for i in range(i):
            loop_vars.append(b[0])
            loop_vars.append(b[0])
            b.remove(b[0])
        loop_to_var.append(loop_vars)
    # list of local polynomials assigned to vertices
    vert_to_Poly = []
    for v in graph.vertices():
        plug_in = []
        for v_edge in graph.edges_incident(v):
            v_edge = tuple( sorted(list(v_edge[:2])) + [v_edge[2]] )
            for var in edge_to_var[v_edge]:
                plug_in.append(var)
        for var in loop_to_var[v]:
            plug_in.append(var)
        used_vars = used_vars + plug_in
        plug_in = plug_in + variables[len(plug_in):]
        Nlocal_at_v = Nlocal(stg.genera[v],valency[v],kappa[v])(*plug_in)
        vert_to_Poly.append(Nlocal_at_v)
    twists = prod(list(set(used_vars)))
    return c*twists*prod(local_poly for local_poly in vert_to_Poly)

def completed_volume(stratum, with_pi=True, verbose=False, one_vertex=False):
    '''
    Return the completed volume of the stratum.
    
    INPUT:
    
    - ``stratum``    -- list, a list of orders of zeros of the stratum
    - ``with_pi``    -- boolean (default `True`), when False returns completed volume as a rational number (volume divided by 
      an appropriate degree of pi
    - ``verbose``    -- boolean (default `False`), when True prints progress of the computation: time to generate and the number       of stable graphs in each codimension and in total; progress of computing contribution of stable graphs
    - ``one_vertex`` -- boolean (default `False`), when True only computed contribution of one-vertex labeled stable graphs.
    
    EXAMPLES:

    Here we compute completed volume of an empty stratum Q(3,1)::
        
        sage: from cvolume import completed_volume 
        sage: completed_volume([3,1])
        23/90*pi^4
        
    Here we demonstrate the verbose mode by computing completed volume of stratum Q(1,-1)::
        
        sage: completed_volume([3, 1, 1, -1], verbose = True)
        Computing completed volume of stratum [3, 1, 1, -1]...
        Generated 2 codimension 1 graphs in ... s
        Generated 4 codimension 2 graphs in ... s
        Generated 3 codimension 3 graphs in ... s
        The total number of stable graphs for stratum [3, 1, 1, -1] is: 9.
        Generated all stable graphs for stratum [3, 1, 1, -1] in: ... s
        Computed contribution of 9/9 graphs. Time elapsed: ... s
        Completed volume of [3, 1, 1, -1] is computed in: ... s
        Completed volume of [3, 1, 1, -1] is: 7/60*pi^6
        7/60*pi^6
        
    Here we compute one-vertex graphs contribution to the completed volume of Q(3,1,1,-1)::
    
        sage: completed_volume([3,1,1,-1], one_vertex=True)
        1346/14175*pi^6
        
    Here are some examples for principal strata, where completed volume coincides with Masur-Veech volume::
    
        sage: assert completed_volume([1,-1,-1,-1,-1,-1]) == 1*pi^4 
        sage: assert completed_volume([1,1,-1,-1]) == 1/3*pi^4
        sage: assert completed_volume([1,-1]) == 2/3*pi^2
        sage: assert completed_volume([-1,-1,-1,-1]) == 2*pi^2
    '''
    def max_weight(stratum): return sum(stratum)/ZZ(2) + len(stratum)/ZZ(2) + stratum.count(-1) + 1
    higher_part = [i for i in stratum if i > 1]
    lower_part = [i for i in stratum if i not in higher_part]
    for new_higher_part in reversed(Combinations(higher_part).list()):
        rest_is_odd = (sum(higher_part) - sum(new_higher_part)) % 2
        w = max_weight(new_higher_part + lower_part) - rest_is_odd
        s_part = tuple(sorted((i+1)/ZZ(2) for i in new_higher_part))
        _ = Fs(s_part,w)
    d = c_d(stratum)
    f = c_f(stratum)
    mu = prod([factorial(stratum.count(i)) for i in range(-1,max(stratum)+1)])        
    if verbose:
        tic_total = time.time()
        print(f"Computing completed volume of stratum {stratum}...")
        tic = time.time()
    stgs = stable_lab_graphs(stratum, one_vertex=one_vertex, verbose=verbose)
    total_num = len(stgs)
    vol = 0
    count, period = 0, 10
    tic = time.time()
    for stg in stgs:      
        vol += operator(graph_poly(stg)) 
        count += 1        
        if verbose and (count%period == 0 or count == total_num):
            toc = time.time()
            print(f"\rComputed contribution of {count}/{total_num} graphs. Time elapsed: {float2time(toc-tic,5)}", end = "") 
    if verbose:
        toc_total = time.time()
        print(f"\nCompleted volume of {stratum} is computed in: {float2time(toc_total-tic_total,5)}")
        print(f"Completed volume of {stratum} is: {f*mu*vol*pi**d}")
    if with_pi: return f*mu*vol*pi**d
    else: return f*mu*vol

def cvolume_by_graphs(stratum, graphs):
    '''
    Return the contribution of given stable graphs to the completed volume of the stratum.
    
    INPUT:
    
    - ``stratum``    -- list, a list of orders of zeros of the stratum
    - ``graphs``     -- iterable, e.g. a set, of Labeled Stable Graphs
    '''
    f = c_f(stratum)
    mu = prod(factorial(stratum.count(i)) for i in range(-1,max(stratum)+1))
    return f*mu*sum(operator(graph_poly(stg) for stg in graphs))

def cvolume_by_cylinders(stratum):    
    '''
    Return the contributions of k-cylinder surfaces to the completed volume of the stratum for each k.
    '''
    stgs = stable_lab_graphs(stratum, by_codim = True)
    for i in range(1,len(stgs)):    # don't include the 0-th set, as it consists of the original non-degenrated graph
        print(f"Volume contribution of {i+1}-cylinder surfaces is {cvolume_by_graphs(stratum, stgs[i])}")
        
def CKazarian(g,k):
    '''
    Return Kazarian constant for recursion on volumes of principal strata.
    '''
    if g<=0 or k<0 or k>g: return 0
    if [g,k] == [1,0]: return ZZ(1)/12
    return ZZ(g-k+1)/(5*g-k-2)*CKazarian(g,k-1) + ZZ(5*g-6-k)*(5*g-4-k)/ZZ(12)*CKazarian(g-1,k)+\
    ZZ(1)/2*sum(CKazarian(g1,k1)*CKazarian(g-g1,k-k1) for g1 in range(1,g) for k1 in range(k+1))

def principal_volume(g_or_stratum,n=-1):
    '''
    Return the Masur-Veech volume of the principal stratum. Input can be either a single argument that is a list of orders of zeroes or two arguments g,n, where g is the genus and n is the number of simple poles. 
    '''
    if g_or_stratum in ZZ:
        assert n != -1, "Input Error: the input must be either g,n or stratum"
        g = g_or_stratum        
    elif type(g_or_stratum) == list:
        assert n == -1, "Input Error: the input must be either g,n or stratum"
        stratum = g_or_stratum
        n = stratum.count(-1)
        g = (sum(stratum)+4)/ZZ(4)
    else: return "Input Error: the input must be either g,n or stratum"       
    if g == 0:
        if n<3: return 0
        else: return pi**(2*n-6)/ZZ(2**(n-5))
    elif [g,n] == [1,0]: return 0
    else: return ZZ(2**(2*g+1))*pi**(6*g-6+2*n)*factorial(4*g-4+n)/factorial(6*g-7+2*n)*sum(dfactorial(5*g-7+2*n-k)/dfactorial(5*g-3-k)*CKazarian(g,k) for k in range(g+1)) 

def asymptotic_volume(stratum):
    '''
    Return conjectural asymptotics for the volume of the stratum from [ADGZZ]. Note that formula is only conjectured in the case when the number of poles is logarithmically small compared to the genus.
    '''
    return ZZ(4)/pi*prod(2**(d+2)/(d+2) for d in stratum)
