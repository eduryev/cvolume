from sage.all import ZZ, diff, factorial, Partitions, Partition, prod
import time
import sys
from .utils import *
from admcycles import psiclass

# t0,t1,t2,t3,t4,t5,t6 = R.gens()[:7]
S = R['z']
z = S('z')
t0,t1,t2,t3,t4,t5,t6 = S('t0*z'),S('t1*z^2'),S('t2*z^3'),S('t3*z^4'),S('t4*z^5'),S('t5*z^6'),S('t6*z^7')
t7,t8,t9,t10,t11,t12,t13 = S('t7*z^8'),S('t8*z^9'),S('t9*z^10'),S('t10*z^11'),S('t11*z^12'),S('t12*z^13'),S('t13*z^14')
t_str = {t0:'t0',t1:'t1',t2:'t2',t3:'t3',t4:'t4',t5:'t5',t6:'t6',t7:'t7',t8:'t8',t9:'t9',t10:'t10',t11:'t11',t12:'t12',t13:'t13'}

def monom(par):
    '''
    Given a partition ``par`` = :math:`\\left[0^{i_0},1^{i_1},\\ldots,n^{i_n}\\right]` return a monomial (in Multivariate Polynomial Ring over Q) :math:`{t_0}^{i_0}\\cdot\\ldots\\cdot{t_n}^{i_n} \\big/ {i_0}!\\cdot\\ldots\\cdot{i_n}!`.
    
    EXAMPLE:
    
    Here we generate all monomials of weight 5::
    
        sage: from cvolume.series import monom
        sage: [monom(l) for l in Partitions(5)]
        [t4, t0*t3, t1*t2, 1/2*t0^2*t2, 1/2*t0*t1^2, 1/6*t0^3*t1, 1/120*t0^5]
    '''
    exp = par.to_exp()
    t_vars = [t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13]
    #print(par)
    return prod(ZZ(1)/factorial(k) for k in exp)*prod([t_vars[i-1] for i in par])

def coeff(par):
    '''
    Return the coefficient of the monomial correspoding to the partition par in the partition function F.
    
    EXAMPLE:
    
    Here are some examples of the coefficients of the partition function in low weight::
    
        sage: from cvolume.series import coeff
        sage: [coeff(Partition([1,1,1])),coeff(Partition([2,2,2])),coeff(Partition([4,4,1])),coeff(Partition([1,1,1,1]))]
        [1, 1/12, 29/2880, 0]
    '''
    if par == Partition([1,1,1]):
        return 1
    elif ((sum([i-1 for i in par])-len(par)+3)/ZZ(3)) in ZZ and ((sum([i-1 for i in par])-len(par)+3)/ZZ(3))>=0:
        n = len(par)
        g = ((sum([i-1 for i in par])-len(par)+3)/ZZ(3)).ceil()
        psi = prod(psiclass(i+1,g,n)**(par[i]-1) for i in range(0,len(par)))
        return psi.evaluate()
    else:
        return ZZ(0)

def diff(F,*args):
    P = F(z=1)
    str2P_var = {str(v):v for v in P.variables()}
    d = {}
    for t in args:
        if t_str[t] in str2P_var:
            d[t] = str2P_var[t_str[t]] 
        else:
            return S.zero()
    P = P.derivative(*(d[t] for t in args))
    return P(t0=t0,t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,t7=t7,t8=t8,t9=t9,t10=t10,t11=t11,t12=t12,t13=t13)

def get_Fs2(F,w):
    return 12*diff(F,t2) - diff(F,t0,t0)/2 - diff(F,t0)._mul_trunc_(diff(F,t0),w)/2

# def get_Fs2(F,w):
#     return 12*diff(F,t2).truncate(w) - diff(F,t0,2).truncate(w)/2 - diff(F,t0).truncate(w)._mul_trunc_(diff(F,t0).truncate(w))/2

def get_Zs2(Z):
    return 12*diff(Z,t2) - diff(Z,t0,2)/2

def get_Fs3(F):
    return 120*diff(F,t3) - 6*diff(F,t0,t1) - 6*diff(F,t1)*diff(F,t0) + 5*diff(F,t0)/4

def get_Zs3(Z):
    return 120*diff(Z,t3) - 6*diff(Z,t0,t1) + 5*diff(Z,t0)/4

def get_Fs4(F):
    return 1680*diff(F, t4) - 18*diff(F, t1, t1) - 18*diff(F, t1)**2 - 60*diff(F, t0, t2) - 60*diff(F, t2)*diff(F, t0) + 7*diff(F, t0, t0, t0)/6 + 7*diff(F, t0)*diff(F, t0, t0)/2 + 7*diff(F, t0)**3/6 + 49*diff(F, t1)/2 - ZZ(35)/96

def get_Fs22(F):
    return ZZ(1)/2*(144*diff(F,t2,t2) - 840*diff(F,t3) - 12*diff(F,t0,t0,t2) - 24*diff(F,t0)*diff(F,t0,t2) + 24*diff(F,t0,t1) + 24*diff(F,t1)*diff(F,t0) + diff(F,t0,4)/4 + diff(F,t0)*diff(F,t0,3) + diff(F,t0,t0)**2/2 + diff(F,t0)**2*diff(F,t0,t0) - 3*diff(F,t0))

def get_Zs22(Z):
    return ZZ(1)/2*(144*diff(Z,t2,t2) - 840*diff(Z,t3) - 12*diff(Z,t0,t0,t2) + 24*diff(Z,t0,t1) + diff(Z,t0,4)/4 - 3*diff(Z,t0))

def get_Fs5(F):
    return 30240*diff(F, t5) - 360*diff(F, t1, t2) - 360*diff(F, t2)*diff(F, t1) - 840*diff(F, t0, t3) - 840*diff(F, t0)*diff(F, t3) + 27*diff(F, t0, t0, t1) + 27*diff(F, t1)*diff(F, t0, t0) + 54*diff(F, t0)*diff(F, t0, t1) + 27*diff(F, t1)*diff(F, t0)*diff(F, t0) + 585*diff(F, t2) - 105*diff(F, t0, t0)/8 - 105*diff(F, t0)*diff(F, t0)/8

def get_Fs23(F):
    return 1440*diff(F, t2, t3) - 15120*diff(F, t4) - 60*diff(F, t0, t0, t3) - 120*diff(F, t0)*diff(F, t0, t3) - 72*diff(F, t0, t1, t2) - 72*diff(F, t1)*diff(F, t0, t2) - 72*diff(F, t1, t2)*diff(F, t0) + 90*diff(F, t1, t1) + 90*diff(F, t1)**2 + 375*diff(F, t0, t2) + 360*diff(F, t2)*diff(F, t0) + 3*diff(F, t0, t0, t0, t1) + 6*diff(F, t0, t1)*diff(F, t0, t0) + 9*diff(F, t0)*diff(F, t0, t0, t1) + 3*diff(F, t1)*diff(F, t0, t0, t0) + 6*diff(F, t1)*diff(F, t0)*diff(F, t0, t0) + 6*diff(F, t0)**2*diff(F, t0, t1) - 45*diff(F, t0, t0, t0)/8 - 65*diff(F, t0)*diff(F, t0, t0)/4 - 5*diff(F, t0)**3 - 165*diff(F, t1)/2 + ZZ(29)/32

def get_Fs6(F):
    return 665280*diff(F, t6) - 1800*diff(F, t2, t2) - 1800*diff(F, t2)**2 - 5040*diff(F, t1, t3) - 5040*diff(F, t3)*diff(F, t1) - 15120*diff(F, t0, t4) - 15120*diff(F, t4)*diff(F, t0) + 16170*diff(F, t3) + 198*diff(F, t0, t1, t1) + 396*diff(F, t1)*diff(F, t0, t1) + 198*diff(F, t1, t1)*diff(F, t0) + 198*diff(F, t1)**2*diff(F, t0) + 330*diff(F, t0, t0, t2) + 330*diff(F, t2)*diff(F, t0, t0) + 660*diff(F, t0)*diff(F, t0, t2) + 330*diff(F, t2)*diff(F, t0)**2 - 33*diff(F, t0, t0, t0, t0)/8 - 33*diff(F, t0)*diff(F, t0, t0, t0)/2 - 99*diff(F, t0, t0)**2/8 - 99*diff(F, t0)**2*diff(F, t0, t0)/4 - 33*diff(F, t0)**4/8 - 891*diff(F, t0, t1)/2 - 891*diff(F, t1)*diff(F, t0)/2 + 1155*diff(F, t0)/32

def get_Fs24(F):
    return 385*diff(F, t0)*diff(F, t0)/8 + 385*diff(F, t0, t0)/8 - 21*diff(F, t0)*diff(F, t0)*diff(F, t0, t0, t0)/4 - 7*diff(F, t0)*diff(F, t0, t0)*diff(F, t0, t0) - 7*diff(F, t0, t0, t0, t0, t0)/12 - 147*diff(F, t0)*diff(F, t0)*diff(F, t1) - 637*diff(F, t0, t0, t1)/4 + 102*diff(F, t0)*diff(F, t0)*diff(F, t0, t2) + 18*diff(F, t0, t1)*diff(F, t0, t1) - 7*diff(F, t0)*diff(F, t0)*diff(F, t0)*diff(F, t0, t0)/2 + 60*diff(F, t2)*diff(F, t0)*diff(F, t0, t0) - 35*diff(F, t0)*diff(F, t0, t0, t0, t0)/12 - 21*diff(F, t0, t0)*diff(F, t0, t0, t0)/4 - 637*diff(F, t0)*diff(F, t0, t1)/2 + 36*diff(F, t1)*diff(F, t0)*diff(F, t0, t1) + 132*diff(F, t0)*diff(F, t0, t0, t2) + 30*diff(F, t2)*diff(F, t0, t0, t0) + 102*diff(F, t0, t2)*diff(F, t0, t0) + 18*diff(F, t0)*diff(F, t0, t1, t1) + 18*diff(F, t1)*diff(F, t0, t0, t1) - 720*diff(F, t0, t2)*diff(F, t2) - 1680*diff(F, t0)*diff(F, t0, t4) - 720*diff(F, t2, t2)*diff(F, t0) - 432*diff(F, t1, t2)*diff(F, t1) - 720*diff(F, t0, t2, t2) - 840*diff(F, t0, t0, t4) + 44*diff(F, t2, t0, t0, t0) + 20160*diff(F, t2, t4) - 216*diff(F, t1, t1, t2) - 147*diff(F, t1)*diff(F, t0, t0) + 9*diff(F, t1, t1, t0, t0) + 2520*diff(F, t2)*diff(F, t1) + 6720*diff(F, t0)*diff(F, t3) + 6720*diff(F, t0, t3) - 2835*diff(F, t2) + 2814*diff(F, t1, t2) - 332640*diff(F, t5)

def get_Fs33(F):
    return ZZ(1)/2*(14400*diff(F, t3, t3) - 332640*diff(F, t5) - 1440*diff(F, t0, t1, t3) - 1440*diff(F, t1, t3)*diff(F, t0) - 1440*diff(F, t0, t3)*diff(F, t1) + 7020*diff(F, t0, t3) + 6720*diff(F, t0)*diff(F, t3) + 2160*diff(F, t1, t2) + 2160*diff(F, t2)*diff(F, t1) + 36*diff(F, t1, t1, t0, t0) + 36*diff(F, t1, t0)*diff(F, t1, t0) + 36*diff(F, t1, t1)*diff(F, t0, t0) + 72*diff(F, t1)*diff(F, t0, t0, t1) + 72*diff(F, t0)*diff(F, t0, t1, t1) + 36*diff(F, t1)*diff(F, t1)*diff(F, t0, t0) + 72*diff(F, t1)*diff(F, t0)*diff(F, t0, t1) + 36*diff(F, t1, t1)*diff(F, t0)*diff(F, t0) - 2400*diff(F, t2) - 165*diff(F, t0, t0, t1) - 165*diff(F, t1)*diff(F, t0, t0) - 150*diff(F, t1)*diff(F, t0)*diff(F, t0) - 315*diff(F, t0)*diff(F, t0, t1) + 725*diff(F, t0, t0)/16 + 175*diff(F, t0)*diff(F, t0)/4)

def get_Fs7(F):
    return 17297280*diff(F,t7) - 50400*diff(F,t2,t3) - 50400*diff(F,t3)*diff(F,t2) - 90720*diff(F,t1,t4) - 90720*diff(F,t4)*diff(F,t1) - 332640*diff(F,t0,t5) - 332640*diff(F,t5)*diff(F,t0) + 468*diff(F,t1,t1,t1) + 1404*diff(F,t1)*diff(F,t1,t1) + 468*diff(F,t1)**3 + 4680*diff(F,t0,t1,t2) + 4680*diff(F,t1)*diff(F,t0,t2) + 4680*diff(F,t2)*diff(F,t0,t1) + 4680*diff(F,t1,t2)*diff(F,t0) + 4680*diff(F,t2)*diff(F,t1)*diff(F,t0) + 5460*diff(F,t0,t0,t3) + 5460*diff(F,t3)*diff(F,t0,t0) + 10920*diff(F,t0)*diff(F,t0,t3) + 5460*diff(F,t3)*diff(F,t0)**2 + 507780*diff(F,t4) - 10725*diff(F,t0,t2) - 10725*diff(F,t2)*diff(F,t0) - 5577*diff(F,t1,t1)/2 - 5577*diff(F,t1)**2/2 - 143*diff(F,t0,t0,t0,t1) - 143*diff(F,t1)*diff(F,t0,t0,t0) - 429*diff(F,t0,t1)*diff(F,t0,t0) - 429*diff(F,t0)*diff(F,t0,t0,t1) - 429*diff(F,t0)**2*diff(F,t0,t1) - 429*diff(F,t1)*diff(F,t0)*diff(F,t0,t0) - 143*diff(F,t1)*diff(F,t0)**3 + 1001*diff(F,t0,t0,t0)/8 + 3003*diff(F,t0)*diff(F,t0,t0)/8 + 1001*diff(F,t0)**3/8 + 27027*diff(F,t1)/16 - ZZ(5005)/384

def get_Fs222(F):
    return ZZ(1)/6*(1728*diff(F,t2,t2,t2) - 30240*diff(F,t2,t3) - 216*diff(F,t0,t0,t2,t2) - 432*diff(F,t0,t2)**2 - 432*diff(F,t0)*diff(F,t0,t2,t2) + 864*diff(F,t0,t1,t2) + 864*diff(F,t1)*diff(F,t0,t2) + 864*diff(F,t1,t2)*diff(F,t0) + 1260*diff(F,t0,t0,t3) + 2520*diff(F,t0)*diff(F,t0,t3) + 9*diff(F,t0,t0,t0,t0,t2) + 36*diff(F,t0,t2)*diff(F,t0,t0,t0) + 36*diff(F,t0)*diff(F,t0,t0,t0,t2) + 36*diff(F,t0,t0)*diff(F,t0,t0,t2) + 72*diff(F,t0)*diff(F,t0,t2)*diff(F,t0,t0) + 36*diff(F,t0)**2*diff(F,t0,t0,t2) + 151200*diff(F,t4) - 576*diff(F,t1,t1) - 576*diff(F,t1)**2 - 2628*diff(F,t0,t2) - 2520*diff(F,t2)*diff(F,t0) - 36*diff(F,t0,t0,t0,t1) - 108*diff(F,t0)*diff(F,t0,t0,t1) - 36*diff(F,t1)*diff(F,t0,t0,t0) - 72*diff(F,t0,t1)*diff(F,t0,t0) - 72*diff(F,t1)*diff(F,t0)*diff(F,t0,t0) - 72*diff(F,t0)**2*diff(F,t0,t1) - (ZZ(1)/8)*diff(F,t0,t0,t0,t0,t0,t0) - (ZZ(3)/4)*diff(F,t0)*diff(F,t0,t0,t0,t0,t0) - (ZZ(3)/2)*diff(F,t0,t0)*diff(F,t0,t0,t0,t0) - (ZZ(5)/4)*diff(F,t0,t0,t0)**2 - 6*diff(F,t0)*diff(F,t0,t0)*diff(F,t0,t0,t0) - diff(F,t0,t0)**3 - (ZZ(3)/2)*diff(F,t0)**2*diff(F,t0,t0,t0,t0) - 3*diff(F,t0)**2*diff(F,t0,t0)**2 - diff(F,t0)**3*diff(F,t0,t0,t0) + (ZZ(63)/2)*diff(F,t0,t0,t0) + 90*diff(F,t0)*diff(F,t0,t0) + 27*diff(F,t0)**3 + 378*diff(F,t1) - ZZ(63)/20)

class PartitionFunctions:
    AC_formulae = {(2,):get_Fs2, (3,):get_Fs3, (4,):get_Fs4, (2,2):get_Fs22, (5,):get_Fs5, (2,3):get_Fs23, (6,):get_Fs6,\
                   (2,4):get_Fs24, (3,3):get_Fs33, (7,):get_Fs7, (2,2,2):get_Fs222}
    Z_formulae = {(2,):get_Zs2}
    
    shifts = {(2,):3, (3,):4, (4,):5, (2,2):6, (5,):6, (2,3):7, (6,):7, (2,4):8, (3,3):8, (7,):8, (2,2,2):9}
    times = {():time_for_F, (2,):time_for_Fs2, (3,):time_for_Fs3, (4,):time_for_Fs4, (2,2):time_for_Fs22, (5,):time_for_Fs5,\
             (2,3):time_for_Fs23, (6,):time_for_Fs6,(2,4):time_for_Fs24, (3,3):time_for_Fs33, (7,):time_for_Fs7,\
             (2,2,2):time_for_Fs222}
    def __init__(self):
        self.F_weights = {}
        self.Z_weights = {}
        self.F_series = {}
        self.Z_series = {}
        self.verbose = False
        
#     def log_partition_function(self,w):
#         Z_max_weight = self.Z_weights.get((),-1)
#         Z = self.Z_series.get((),R.zero())
#         if w > Z_max_weight:
#             F = 
#             Z = F.log(prec = w)
#             self.Z_weights[()] = w
#             self.Z_series[()] = Z
#             toc = time.time()
#             if self.verbose: print(f"    Done updating the partition function F from max_weight {F_max_weight} to {w} in: {float2time(toc-tic,2)}")
#         return F
    
    def partition_function(self,w):
        '''
        Return the partition function F truncated at weight w.
        
        EXAMPLE:
        
        The partition function truncated at weight 10::
        
            sage: from cvolume import Fs
            sage: Fs.partition_function(10)
            1/6*t0^3*t1^3 + 1/8*t0^4*t1*t2 + 1/120*t0^5*t3 + 1/6*t0^3*t1^2 + 1/120*t1^5 + 1/24*t0^4*t2 + 1/6*t0*t1^3*t2 + 1/6*t0^2*t1*t2^2 + 1/8*t0^2*t1^2*t3 + 7/144*t0^3*t2*t3 + 1/36*t0^3*t1*t4 + 1/576*t0^4*t5 + 1/6*t0^3*t1 + 1/96*t1^4 + 1/8*t0*t1^2*t2 + 1/24*t0^2*t2^2 + 1/16*t0^2*t1*t3 + 1/144*t0^3*t4 + 1/6*t0^3 + 1/72*t1^3 + 1/12*t0*t1*t2 + 7/1440*t2^3 + 1/48*t0^2*t3 + 29/1440*t1*t2*t3 + 29/5760*t0*t3^2 + 1/192*t1^2*t4 + 11/1440*t0*t2*t4 + 1/288*t0*t1*t5 + 1/2304*t0^2*t6 + 1/48*t1^2 + 1/24*t0*t2 + 29/5760*t2*t3 + 1/384*t1*t4 + 607/2903040*t4^2 + 1/1152*t0*t5 + 503/1451520*t3*t5 + 77/414720*t2*t6 + 5/82944*t1*t7 + 1/82944*t0*t8 + 1/24*t1 + 1/1152*t4 + 1/82944*t7
        '''
        F_max_weight = self.F_weights.get((),-1)
        F = self.F_series.get((),R.zero()) 
        if w > F_max_weight:
            tic = time.time()
            time_est = time_for_F(w) - time_for_F(F_max_weight)
            if time_est > 120:
                command = input(f"    Partition function update might take more than {float2time(time_est,2)}. Do you want to \
                continue? Print 'n' to abort, to continue press Enter.")
                if command == 'n':
                    sys.exit("The computation was aborted by user due to time constraints.")
            if self.verbose: print(f"    Updating the partition function F from max_weight {F_max_weight} to {w}...Estimated time: {float2time(time_est,2)}")
            for i in range(max(F_max_weight,0)+1,w+1):
                F += sum(coeff(par)*monom(par) for par in Partitions(i))
            self.F_weights[()] = w
            self.F_series[()] = F
            toc = time.time()
            if self.verbose: print(f"    Finished in: {float2time(toc-tic,2)}")
        return F
    
    def __call__(self, w, s_part=None):
        if s_part is None:
            return self.partition_function(w)
        Fs_max_weight = self.F_weights.get(s_part,-1)
        Fs = self.F_series.get(s_part,R.zero())             
        if w > Fs_max_weight:
            tic = time.time()
            time_est = self.times[s_part](w)
            if time_est > 120:
                command = input(f"Fs function update for s = {s_part}, w = {w} might take more than {float2time(time_est,2)}.\
                Do you want to continue? Print 'n' to abort, to continue press Enter..")              
                if command == 'n':
                    sys.exit("The computation was aborted by user due to time constraints.")
            if self.verbose: print(f"Updating Fs function for s = {s_part} from max_weight {Fs_max_weight} to {w}...Estimated time: {float2time(time_est,2)}")
            F = self.partition_function(w+self.shifts[s_part])
            Fs = self.AC_formulae[s_part](F,w+1)
            self.F_weights[s_part] = w
            self.F_series[s_part] = Fs
            toc = time.time()
            if self.verbose: print(f"Finished in: {float2time(toc-tic,2)}")     
        return Fs
    
    def reset(self):
        self.F_weights = {}
        self.Z_weights = {}
        self.F_series = {}
        self.Z_series = {}

    
Fs = PartitionFunctions()
