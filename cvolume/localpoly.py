from sage.all import ZZ, factorial, prod, Partitions, Permutations, OrderedSetPartitions, diff
from .utils import *
from .series import Fs
        
def stratum_to_F(g,n,stratum):
    '''
    Return partition function Fs corresponding to the stratum. Function is truncated at a minimal admissible weight.
    '''
    s_part = tuple(sorted((i+1)/ZZ(2) for i in stratum if i > 1))
    w = sum(stratum)/ZZ(4) + len(stratum)/ZZ(2) + stratum.count(-1) + n/ZZ(2)
    return Fs(w,s_part)

def vanish(variables,indices):
    '''
    Return a copy of the list of variables, in which elements at given indices are replaces by 0.
    '''
    l = variables[:]
    for i in indices:
        l[i-1] = 0
    return l

def shift(variables,k):
    '''
    Return a copy of the list of variables shifted k units to the left.
    '''
    l = variables[:]
    for i in range(len(l)-k):
        l[i] = l[i+k]
    for i in range(len(l)-k,len(l)):
        l[i] = 0
    return l

def Nlocal(g,n,stratum,labeled=False,mode='derivative'):
    '''
    Return a local polynomial :math:`N_{g,n}^{\\kappa}(b_1, ..., b_n)`, where :math:`\\kappa` is stratum.
    
    INPUT:
    
    - ``g``       -- int, genus
    - ``n``       -- int, number of boundaries
    - ``stratum`` -- list or tuple, orders of zeros of the stratum
    - ``labeled`` -- boolean (default `False`), whether to label zeros or not
    - ``mode``    -- 'derivative' (default) or 'recursive', 'derivative' uses Arbarello-Cornalba formulae and 'recursive' uses recursion on local polynomials dervied from the same AC formulae
    
    OUTPUT:
    
    - a symmetric Sage Multivariate Polynomial in variables b1, b2, ... .
    
    EXAMPLES:
    
    Here we compute local polynomial of genus 2, two bondary components and stratum :math:`\\mathcal{Q}(5, 1, 1, 1)`::
    
        sage: from cvolume import Nlocal
        sage: S = PolynomialRing(QQ,['b%d' % i for i in range(1,10)])
        sage: b1,b2,b3,b4,b5 = S.gens()[:5]
        sage: Nlocal(2, 2, (5, 1, 1, 1))
        1/192*b1^6 + 95/3072*b1^4*b2^2 + 95/3072*b1^2*b2^4 + 1/192*b2^6
   
    And here is the same local polynomial, but with labeled zeros::
    
        sage: Nlocal(2, 2, (5, 1, 1, 1), labeled = True)
        1/32*b1^6 + 95/512*b1^4*b2^2 + 95/512*b1^2*b2^4 + 1/32*b2^6
    
    Another example for a different stratum::
    
        sage: Nlocal(2, 2, (7, 1))
        35/192*b1^4 + 35/64*b1^2*b2^2 + 35/192*b2^4
    
    Here we show that `recursive` and `derivative` methods give consistent answers::
    
        sage: assert Nlocal(0, 5, [3, 1, 1, 1], mode = 'recursive') == Nlocal(0, 5, [3, 1, 1, 1])
        sage: assert Nlocal(1, 3, [3, 1, 1, 1], mode = 'recursive') == Nlocal(1, 3, [3, 1, 1, 1])
        sage: assert Nlocal(0, 4, [3, 1, 1, -1], mode = 'recursive') == Nlocal(0, 4, [3, 1, 1, -1])
        sage: assert Nlocal(0, 3, [3, 1, -1, -1], mode = 'recursive') == Nlocal(0, 3, [3, 1, -1, -1])
    '''
    if not stratum:
        stratum = [1]*(4*g-4+2*n)
    cache_poly = {}
    cache_labeling = {}
    stratum = tuple(stratum)
    b = list(S.gens())
    def memo_Nlocal(g,n,stratum,labeled,mode):
        #print('Computing for (%s,%s,%s)...Cache size %s' % (g,n,stratum,len(cache_poly)))
        if not stratum or n != 2-2*g+1/ZZ(2)*sum(stratum) or g<0 or n<1: return S.zero()
        if type(stratum) == list: stratum = tuple(stratum)
        if (g,n,stratum) in cache_poly:
            if labeled: return cache_poly[(g,n,stratum)]
            else: return cache_poly[(g,n,stratum)]/cache_labeling[(g,n,stratum)] 
        if not labeled:
            memo_Nlocal(g,n,stratum,True,mode)
            return cache_poly[(g,n,stratum)]/cache_labeling[(g,n,stratum)]
        #weight_check(g,n,stratum)
        m = [stratum.count(2*i-1) for i in range((max(stratum)+1)/ZZ(2)+1)]
        if -1 in stratum:                    # case with poles
            elderstratum = list(stratum)
            elderstratum.remove(-1)
            elderstratum.append(1)
            N_lab = memo_Nlocal(g,n+1,elderstratum,True,mode)(*vanish(b,[n+1]))
            for i in range(len(elderstratum)):
                if elderstratum[i]>1:
                    newstratum = list(elderstratum)
                    newstratum[i] = newstratum[i]-2
                    N_lab -= elderstratum[i]*memo_Nlocal(g,n,newstratum,True,mode)
            N_lab *= ZZ(1)/elderstratum.count(1)
            cache_poly[(g,n,stratum)] = N_lab
            cache_labeling[(g,n,stratum)] = ZZ(prod(factorial(stratum.count(i)) for i in range(-1,max(stratum)+1)))
            return cache_poly[(g,n,stratum)]
        elif mode == 'derivative' or len(m) == 2:             # case without poles
            _Fs = stratum_to_F(g,n,stratum) 
            deg = ZZ(3*g-3+n-1/2*sum(d-1 for d in stratum))
            M = sum(m[i]*(i-1) for i in range(len(m)))
            mondeg = Partitions(deg+n, length=n)
            const = 2**(5*g-6+2*n-2*M)
            labeling = prod(factorial(stratum.count(i)) for i in range(1,max(stratum)+1))
            N_lab = ZZ(1)/const*labeling*sum(prod(ZZ(1)/factorial(d-1) for d in par)*\
                            prod(factorial(i) for i in par.to_exp())*\
                            _Fs.monomial_coefficient(prod(R.gen(d-1) for d in par))*\
                            sum(prod(S.gen(j)**(2*(sympar[j]-1)) for j in range(n)) for sympar in Permutations(par))\
                            for par in mondeg)
            cache_poly[(g,n,stratum)] = N_lab
            cache_labeling[(g,n,stratum)] = ZZ(prod(factorial(stratum.count(i)) for i in range(-1,max(stratum)+1)))
            return cache_poly[(g,n,stratum)]
        elif mode == 'recursive':
            assert m[2] == 1 and len(m) == 3, "'recursive' mode is only available for stratum = [3,1,1,...,-1,-1]"  
            first_term = 2**4*diff(memo_Nlocal(g,n+1,[1]*(m[1]+5),False,mode), b[n], 4)(*vanish(b,[n+1]))
            second_term = -ZZ(1)/2*2*memo_Nlocal(g-1,n+2,[1]*(m[1]+3),False,mode)(*vanish(b,[n+1,n+2]))
            third_term = 0
            for par in OrderedSetPartitions(b[:n], 2):
                n_1, n_2 = len(par[0]), len(par[1])
                assert n_1 + n_2 == n
                third_term += -ZZ(1)/2*sum(memo_Nlocal(g_1,n_1+1,[1]*(4*g_1-4+2*(n_1+1)),False,mode)(*list(par[0])+[0]*(30-n_1))\
                        *memo_Nlocal(g-g_1,n_2+1,[1]*(4*(g-g_1)-4+2*(n_2+1)),False,mode)(*list(par[1])+[0]*(30-n_2))\
                                  for g_1 in range(g+1))        
            N_unlab = first_term + second_term + third_term            
            cache_labeling[(g,n,stratum)] = ZZ(prod(factorial(stratum.count(i)) for i in range(-1,max(stratum)+1)))
            cache_poly[(g,n,stratum)] = N_unlab*cache_labeling[(g,n,stratum)]
            return cache_poly[(g,n,stratum)]
    memo_Nlocal(g,n,stratum,labeled,mode)
    if labeled: return cache_poly.get((g,n,stratum),S.zero())
    else: return cache_poly.get((g,n,stratum),S.zero())/cache_labeling.get((g,n,stratum),1)
    