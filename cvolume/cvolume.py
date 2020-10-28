from sage.all import factorial, Rational, prod, Combinations, zeta, pi
import itertools
import time
from .localpoly import Nlocal
from .utils import *
from .stable_graphs import *
from .series import Fs

import pickle
with open("cvolume/Elise_volumes.pkl", "rb") as f:
    volumes_dict = pickle.load(f)
    
cvolumes_dict = {}
conjvolumes_dict = {}

def replmon(n):
    '''
    Return the result of application of :math:`\\mathcal{Z}`-operator to the monomial :math:`b_i^n`.
    '''
    if n==0: return 1
    else: return factorial(n)*Rational(zeta(n+1)/pi**(n+1))
    
def operator(Poly):
    '''
    Return the result of application of :math:`\\mathcal{Z}`-operator to the polynomial.
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
    
        sage: stg = LabeledStableGraph([(0, 1, 1)], [1, 1], [[3, -1], [3, -1]])
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
    - ``one_vertex`` -- boolean (default `False`), when True only computes contribution of one-vertex labeled stable graphs.
    
    EXAMPLES:

    Here we compute completed volume of an empty stratum :math:`\\mathcal{Q}(3,1)`::
        
        sage: from cvolume import completed_volume 
        sage: completed_volume([3, 1])
        23/90*pi^4
        
    Here we demonstrate the verbose mode by computing completed volume of stratum :math:`\\mathcal{Q}(1,-1)`::
        
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
        
    Here we compute one-vertex graphs contribution to the completed volume of :math:`\\mathcal{Q}(3,1,1,-1)`::
    
        sage: completed_volume([3, 1, 1, -1], one_vertex=True)
        1346/14175*pi^6
        
    Here are some examples for principal strata, where completed volume coincides with Masur-Veech volume::
    
        sage: assert completed_volume([1, -1, -1, -1, -1, -1]) == 1*pi^4 
        sage: assert completed_volume([1, 1, -1, -1]) == 1/3*pi^4
        sage: assert completed_volume([1, -1]) == 2/3*pi^2
        sage: assert completed_volume([-1, -1, -1, -1]) == 2*pi^2
    '''
    stratum = tuple(sorted(stratum,reverse=True))
    if stratum in cvolumes_dict and one_vertex == False:
        cvolume = cvolumes_dict[stratum]
        if with_pi: return cvolume
        else: return cvolume.coefficient(pi**cvolume.degree(pi))
    def max_weight(stratum): return sum(stratum)/ZZ(2) + len(stratum)/ZZ(2) + stratum.count(-1) + 1
    higher_part = [i for i in stratum if i > 1]
    lower_part = [i for i in stratum if i not in higher_part]
    for new_higher_part in reversed(Combinations(higher_part).list()):
        rest_is_odd = (sum(higher_part) - sum(new_higher_part)) % 2
        w = max_weight(new_higher_part + lower_part) - rest_is_odd
        s_part = tuple(sorted((i+1)/ZZ(2) for i in new_higher_part))
        _ = Fs(w,s_part)
    d = c_d(stratum)
    f = c_f(d)
    mu = prod([factorial(stratum.count(i)) for i in range(-1,max(stratum)+1)])        
    if verbose:
        tic_total = time.time()
        print(f"Computing completed volume of stratum {stratum2print(stratum)}...")
        tic = time.time()
    stgs = stable_lab_graphs(stratum, one_vertex=one_vertex, verbose=verbose)
    total_num = len(stgs)
    vol, count, period = 0, 0, 10
    tic = time.time()
    for stg in stgs:      
        vol += operator(graph_poly(stg)) 
        count += 1        
        if verbose and (count%period == 0 or count == total_num):
            toc = time.time()
            print(f"\rComputed contribution of {count}/{total_num} graphs. Time elapsed: {float2time(toc-tic,3)}. ETA: {float2time((toc-tic)/count*(total_num-count),3)}", end = "") 
    if verbose:
        toc_total = time.time()
        print(f"\nCompleted volume of {stratum2print(stratum)} is: {f*mu*vol*pi**d}. Computed in: {float2time(toc_total-tic_total,3)}")
    cvolume, cvolume_pi = f*mu*vol, f*mu*vol*pi**d
    if one_vertex == False: cvolumes_dict[stratum] = cvolume_pi
    if with_pi: return cvolume_pi
    else: return cvolume

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
        
def CKazarian(g,k,memo=None):
    '''
    Return Kazarian constant for recursion on volumes of principal strata.
    '''
    if g<=0 or k<0 or k>g: return 0
    if (g,k) == (1,0): return ZZ(1)/12
    if (g,k) in memo: return memo[(g,k)]
    const = ZZ(g-k+1)/(5*g-k-2)*CKazarian(g,k-1,memo) + ZZ(5*g-6-k)*(5*g-4-k)/ZZ(12)*CKazarian(g-1,k,memo)+\
    ZZ(1)/2*sum(CKazarian(g1,k1,memo)*CKazarian(g-g1,k-k1,memo) for g1 in range(1,g) for k1 in range(k+1))
    memo[(g,k)] = const
    return const

def principal_volume(g_or_stratum,n=None):
    '''
    Return the Masur-Veech volume of the principal stratum. By convention an empty stratum has volume 0. Input can be either a single argument that is a list of orders of zeroes or two arguments g and n, where g is the genus and n is the number of simple poles. 
    '''
    if g_or_stratum in ZZ:
        assert n != None, "Input Error: the input must be either g,n or stratum"
        g = ZZ(g_or_stratum)
        stratum = [1]*(4*g-4+n) + [-1]*n
    elif type(g_or_stratum) in {list,tuple}:
        assert n == None, "Input Error: the input must be either g,n or stratum"
        stratum = g_or_stratum
        n = stratum.count(-1)
        g = (sum(stratum)+4)/ZZ(4)
    else:
        raise ValueError("Input Error: the input must be either g,n or stratum")
    if sum(stratum) < -4 or sum(stratum)%4 != 0 or n < 0 or (g,n) in {(0,0),(0,1),(0,2),(0,3),(1,0)}: return 0
    memo = {}
    if g == 0: volume = 2*pi**(2*n-6)/ZZ(2**(n-4))
    else: volume = ZZ(2**(2*g+1))*pi**(6*g-6+2*n)*factorial(4*g-4+n)/factorial(6*g-7+2*n)*sum(dfactorial(5*g-7+2*n-k)/dfactorial(5*g-3-k)*CKazarian(g,k,memo) for k in range(g+1)) 
    return volume
    
def MV_volume(stratum,verbose=False,mode='default'):
    '''
    TODO: Add assertions for MV_volume.
    Return the Masur-Veech volume of any startum of quadratic differentials, if it is known. This inlcudes all principal strata, strata covered by the table of volumes of Elise Goujard and strata [3,1,1...-1,-1] and [5,1,1...-1,-1], for which we make use of completed volumes. The mode can be set to 'conjecture', the Masur-Veech volumes are computed using their conjectural relations with completed volumes.
    '''
    stratum = tuple(sorted(stratum,reverse=True))
    if sum(stratum) < -4 or sum(stratum)%4 != 0:
        print(f'Warning: {stratum2print(stratum)} is an empty stratum. Its sum of the orders of zeros is either < -4 or not divisible by 4.')
        return 0
    elif stratum in {(3,1),(1,-1)}:
        return 0
    elif stratum in volumes_dict and mode in {'default','both'}:
        if type(volumes_dict[stratum]) == tuple:
            p,q,deg = volumes_dict[stratum]
            mv_volume = ZZ(p)/q*pi**deg
        else:
            mv_volume = volumes_dict[stratum]
        return mv_volume
    tic = time.time()
    if verbose: print(f"Computing Masur-Veech volume of stratum {stratum2print(stratum)}...")
    d,g,p = c_d(stratum),sum(stratum)/4+1,stratum.count(-1) # dimension, genus, number of poles
    if set(stratum).issubset({1,-1}):
        mv_volume = principal_volume(stratum)  
    elif set(stratum).issubset({3,1,-1}) and stratum.count(3) == 1:        
        boundary = c_f(d)/(c_f(2)*c_f(d-2))*2/3*pi**2*principal_volume(g-1,p+1)
        cvolume = completed_volume(stratum,verbose=verbose)
        mv_volume = cvolume - boundary
        if mode in {'default','both'}: volumes_dict[stratum] = mv_volume
    elif set(stratum).issubset({5,1,-1}) and stratum.count(5) == 1:
        boundary = c_f(d)/(c_f(2)*c_f(d-2))*2/3*pi**2*principal_volume(g-1,p)
        cvolume = completed_volume(stratum,verbose=verbose)
        mv_volume = cvolume - 3*boundary
        if mode in {'default','both'}: volumes_dict[stratum] = mv_volume
    elif mode == 'conjecture' or mode == 'both':
        if stratum in conjvolumes_dict:
            return conjvolumes_dict[stratum]            
        if len(stratum) == 2:
            raise ValueError('There is no conjecture for Masur-Veech volume for this stratum. Stratum needs to have at least three zeros.')
        H0,H2,H4 = ZZ(2)/3*pi**2,ZZ(1)/15*pi**4,ZZ(61)/3402*pi**6
        higher_part = sorted([i for i in stratum if i != 1 and i != -1],reverse=True)
        lower_part = sorted([i for i in stratum if i == 1 or i == -1],reverse=True)
        if higher_part == [7]:
            H0bdry = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [3],verbose=verbose,mode='both')
            H2bdry = c_f(d)/(c_f(4)*c_f(d-4))*H2*MV_volume(lower_part + [-1],verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0*H0*MV_volume(lower_part + [-1],verbose=verbose,mode='both')
            boundary = 5*H0bdry + 3*H2bdry + ZZ(7)/2*H0H0bdry
        elif higher_part == [9]:
            H0bdry = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [5],verbose=verbose,mode='both')
            H2bdry = c_f(d)/(c_f(4)*c_f(d-4))*H2*MV_volume(lower_part + [1],verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0*H0*MV_volume(lower_part + [1],verbose=verbose,mode='both')
            boundary = 7*H0bdry + 9*H2bdry + ZZ(27)/2*H0H0bdry
        elif higher_part == [3,3]:
            H0bdry = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [3,-1],verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0*H0*MV_volume(lower_part + [-1,-1],verbose=verbose,mode='both')
            boundary = 2*H0bdry + H0H0bdry
        elif higher_part == [5,3]:
            H0bdry_5 = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [-1,5],verbose=verbose,mode='both')
            H0bdry_3 = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [1,3],verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0*H0*MV_volume(lower_part + [1,-1],verbose=verbose,mode='both')
            boundary = 3*H0bdry_3 + H0bdry_5 + 3*H0H0bdry
        elif higher_part == [5,5]:
            H0bdry = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [5,1],verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0*H0*MV_volume(lower_part + [1,1],verbose=verbose,mode='both')
            boundary = 6*H0bdry + 9*H0H0bdry        
        elif higher_part == [7,3]:
            H0bdry_33 = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [3,3],verbose=verbose,mode='both')
            H2bdry = c_f(d)/(c_f(4)*c_f(d-4))*H2*MV_volume(lower_part + [3,-1],verbose=verbose,mode='both')
            H0bdry_7 = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume(lower_part + [7,-1],verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0*H0*MV_volume(lower_part + [3,-1],verbose=verbose,mode='both')
            H0H2bdry = c_f(d)/(c_f(2)*c_f(4)*c_f(d-6))*H0*H2*MV_volume(lower_part + [-1,-1],verbose=verbose,mode='both')
            H0H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(2)*c_f(d-6))*H0*H0*H0*MV_volume(lower_part + [-1,-1],verbose=verbose,mode='both')
            boundary = 5*H0bdry_33 + 3*H2bdry + H0bdry_7 + ZZ(17)/2*H0H0bdry + 3*H0H2bdry + ZZ(7)/2*H0H0H0bdry
        elif higher_part == [11]:
            H0bdry = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume([7] + lower_part,verbose=verbose,mode='both')
            H2bdry = c_f(d)/(c_f(4)*c_f(d-4))*H2*MV_volume([3] + lower_part,verbose=verbose,mode='both')
            H4bdry = c_f(d)/(c_f(6)*c_f(d-6))*H4*MV_volume([-1] + lower_part,verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0**2*MV_volume([3] + lower_part,verbose=verbose,mode='both')
            H0H2bdry = c_f(d)/(c_f(2)*c_f(4)*c_f(d-6))*H0*H2*MV_volume([-1] + lower_part,verbose=verbose,mode='both')
            H0H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(2)*c_f(d-6))*H0**3*MV_volume([-1] + lower_part,verbose=verbose,mode='both')
            boundary = 9*H0bdry + 15*H2bdry + 5*H4bdry + 55/2*H0H0bdry + 33*H0H2bdry + 33/2*H0H0H0bdry
        elif higher_part == [3,3,3]:
            H0bdry = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume([3,3,-1] + lower_part,verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0**2*MV_volume([3,-1,-1] + lower_part,verbose=verbose,mode='both')       
            H0H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(2)*c_f(d-6))*H0**3*MV_volume([-1,-1,-1] + lower_part,verbose=verbose,mode='both') 
            boundary = 3*H0bdry + 3*H0H0bdry + H0H0H0bdry
        elif higher_part == [13]:
            print('WARNING: This conjectural formula has not been sufficiently verified.')
            A, B, C = MV_volume([9] + lower_part), MV_volume([5]+lower_part), principal_volume_Kazarian(g-3,poles)
            H0bdry = c_f(d)/(c_f(2)*c_f(d-2))*H0*MV_volume([9]+lower_part,verbose=verbose,mode='both')
            H2bdry = c_f(d)/(c_f(4)*c_f(d-4))*H2*MV_volume([5]+lower_part,verbose=verbose,mode='both')
            H4bdry = c_f(d)/(c_f(6)*c_f(d-6))*H4*MV_volume([1]+lower_part,verbose=verbose,mode='both')
            H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(d-4))*H0**2*MV_volume([5]+lower_part,verbose=verbose,mode='both')
            H0H2bdry = c_f(d)/(c_f(2)*c_f(4)*c_f(d-6))*H0*H2*MV_volume([1]+lower_part,verbose=verbose,mode='both')
            H0H0H0bdry = c_f(d)/(c_f(2)*c_f(2)*c_f(2)*c_f(d-6))*H0**3*MV_volume([1]+lower_part,verbose=verbose,mode='both')
            boundary = 11*H0bdry + 21*H2bdry + 15*H4bdry + 91/2*H0H0bdry + 117*H0H2bdry + 143/2*H0H0H0bdry
        else:
            raise ValueError('There is no conjecture for Masur-Veech volume for this stratum.')
        cvolume = completed_volume(stratum,verbose=verbose)
        mv_volume = cvolume - boundary 
        conjvolumes_dict[stratum] = mv_volume
    else:
        raise ValueError('Masur-Veech volume for this stratum is not known.')    
    if verbose:
        toc = time.time()
        print(f"Masur-Veech volume of {stratum2print(stratum)} is: {mv_volume}. Computed in: {float2time(toc-tic,3)}")
    return mv_volume

def asymptotic_volume(stratum):
    '''
    Return conjectural asymptotics for the volume of the stratum from [ADGZZ]. Note that formula is only conjectured in the case when the number of poles is logarithmically small compared to the genus.
    '''
    return ZZ(4)/pi*prod(2**(d+2)/(d+2) for d in stratum)

