from sage.all import ZZ, diff, factorial, Partitions, Partition, prod
import time
import sys
from .utils import *
from admcycles import psiclass


TZ = T['z']
z = TZ('z')
for i in range(NUM_T_VAR):
    globals()[f't{i}'] = TZ(f't{i}*z^{i+1}')
    globals()[f'T{i}'] = T(f't{i}')
t_vars = [globals()[f't{i}'] for i in range(NUM_T_VAR)]
T_vars = [globals()[f'T{i}'] for i in range(NUM_T_VAR)]
t2T = {globals()[f't{i}']:globals()[f'T{i}'] for i in range(NUM_T_VAR)}
T2t = {globals()[f'T{i}']:globals()[f't{i}'] for i in range(NUM_T_VAR)}

def diff(F,*args):
    '''
    Take polynomial in T-ring and return its derivative in TZ-ring.
    
    '''
    P = F.derivative(*(t2T[t] for t in args))
    if P == 0:
        return TZ.zero()
    else:
        return P(*(t_vars[i] for i in range(NUM_T_VAR)))

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
    return prod(ZZ(1)/factorial(k) for k in exp)*prod([T_vars[i-1] for i in par])

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

def get_Fs2(F,w):
    Fs2 = 12*diff(F,t2).truncate(w+1) - diff(F,t0,t0).truncate(w+1)/2 - diff(F,t0)._mul_trunc_(diff(F,t0),w+1)/2
    return Fs2(z=1)

def get_Zs2(Z):
    return 12*diff(Z,t2) - diff(Z,t0,2)/2

def get_Fs3(F,w):
    Fs3 = 120*diff(F,t3).truncate(w+1) - 6*diff(F,t0,t1).truncate(w+1) - 6*diff(F,t1)._mul_trunc_(diff(F,t0),w+1) + 5*diff(F,t0).truncate(w+1)/4
    return Fs3(z=1)

def get_Zs3(Z):
    return 120*diff(Z,t3) - 6*diff(Z,t0,t1) + 5*diff(Z,t0)/4

def get_Fs4(F,w):
    Fs4 = 1680*diff(F, t4).truncate(w+1) - 18*diff(F, t1, t1).truncate(w+1) - 18*diff(F, t1).power_trunc(2,w+1) - 60*diff(F, t0, t2).truncate(w+1) - 60*diff(F, t2)._mul_trunc_(diff(F, t0),w+1) + 7*diff(F, t0, t0, t0).truncate(w+1)/6 + 7*diff(F, t0)._mul_trunc_(diff(F, t0, t0),w+1)/2 + 7*diff(F, t0).power_trunc(3,w+1)/6 + 49*diff(F, t1).truncate(w+1)/2 - ZZ(35)/96
    return Fs4(z=1)

def get_Fs22(F,w):
    Fs22 = ZZ(1)/2*(144*diff(F,t2,t2).truncate(w+1) - 840*diff(F,t3).truncate(w+1) - 12*diff(F,t0,t0,t2).truncate(w+1) - 24*diff(F,t0)._mul_trunc_(diff(F,t0,t2),w+1) + 24*diff(F,t0,t1).truncate(w+1) + 24*diff(F,t1)._mul_trunc_(diff(F,t0),w+1) + diff(F,t0,t0,t0,t0).truncate(w+1)/4 + diff(F,t0)._mul_trunc_(diff(F,t0,t0,t0),w+1) + diff(F,t0,t0).power_trunc(2,w+1)/2 + diff(F,t0).power_trunc(2,w+1)._mul_trunc_(diff(F,t0,t0),w+1) - 3*diff(F,t0).truncate(w+1))
    return Fs22(z=1)

def get_Zs22(Z):
    return ZZ(1)/2*(144*diff(Z,t2,t2) - 840*diff(Z,t3) - 12*diff(Z,t0,t0,t2) + 24*diff(Z,t0,t1) + diff(Z,t0,4)/4 - 3*diff(Z,t0))

def get_Fs5(F,w):
    Fs5 = 30240*diff(F, t5).truncate(w+1) - 360*diff(F, t1, t2).truncate(w+1) - 360*diff(F, t2)._mul_trunc_(diff(F, t1),w+1) - 840*diff(F, t0, t3).truncate(w+1) - 840*diff(F, t0)._mul_trunc_(diff(F, t3),w+1) + 27*diff(F, t0, t0, t1).truncate(w+1) + 27*diff(F, t1)._mul_trunc_(diff(F, t0, t0),w+1) + 54*diff(F, t0)._mul_trunc_(diff(F, t0, t1),w+1) + 27*diff(F, t1)._mul_trunc_(diff(F, t0).power_trunc(2,w+1),w+1) + 585*diff(F, t2).truncate(w+1) - 105*diff(F, t0, t0).truncate(w+1)/8 - 105*diff(F, t0).power_trunc(2,w+1)/8
    return Fs5(z=1)

def get_Fs23(F,w):
    Fs23 = 1440*diff(F,t2,t3).truncate(w+1) - 15120*diff(F,t4).truncate(w+1) - 60*diff(F,t0,t0,t3).truncate(w+1) - 120*diff(F,t0)._mul_trunc_(diff(F,t0,t3),w+1) - 72*diff(F,t0,t1,t2).truncate(w+1) - 72*diff(F,t1)._mul_trunc_(diff(F,t0,t2),w+1) - 72*diff(F,t1,t2)._mul_trunc_(diff(F,t0),w+1) + 90*diff(F,t1,t1).truncate(w+1) + 90*diff(F,t1).power_trunc(2,w+1) + 375*diff(F,t0,t2).truncate(w+1) + 360*diff(F,t2)._mul_trunc_(diff(F,t0),w+1) + 3*diff(F,t0,t0,t0,t1).truncate(w+1) + 6*diff(F,t0,t1)._mul_trunc_(diff(F,t0,t0),w+1) + 9*diff(F,t0)._mul_trunc_(diff(F,t0,t0,t1),w+1) + 3*diff(F,t1)._mul_trunc_(diff(F,t0,t0,t0),w+1) + 6*diff(F,t1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t0),w+1) + 6*diff(F,t0).power_trunc(2,w+1)._mul_trunc_(diff(F,t0,t1),w+1) - 45*diff(F,t0,t0,t0).truncate(w+1)/8 - 65*diff(F,t0)._mul_trunc_(diff(F,t0,t0),w+1)/4 - 5*diff(F,t0).power_trunc(3,w+1) - 165*diff(F,t1).truncate(w+1)/2 + ZZ(29)/32
    return Fs23(z=1)

def get_Fs6(F,w):
    Fs6 = 665280*diff(F, t6).truncate(w+1) - 1800*diff(F, t2, t2).truncate(w+1) - 1800*diff(F, t2)._mul_trunc_(diff(F, t2),w+1) - 5040*diff(F, t1, t3).truncate(w+1) - 5040*diff(F, t3)._mul_trunc_(diff(F, t1),w+1) - 15120*diff(F, t0, t4).truncate(w+1) - 15120*diff(F, t4)._mul_trunc_(diff(F, t0),w+1) + 16170*diff(F, t3).truncate(w+1) + 198*diff(F, t0, t1, t1).truncate(w+1) + 396*diff(F, t1)._mul_trunc_(diff(F, t0, t1),w+1) + 198*diff(F, t1, t1)._mul_trunc_(diff(F, t0),w+1) + 198*diff(F, t1).power_trunc(2,w+1)._mul_trunc_(diff(F, t0),w+1) + 330*diff(F, t0, t0, t2).truncate(w+1) + 330*diff(F, t2)._mul_trunc_(diff(F, t0, t0),w+1) + 660*diff(F, t0)._mul_trunc_(diff(F, t0, t2),w+1) + 330*diff(F, t2)._mul_trunc_(diff(F, t0),w+1)._mul_trunc_(diff(F, t0),w+1) - 33*diff(F, t0, t0, t0, t0).truncate(w+1)/8 - 33*diff(F, t0)._mul_trunc_(diff(F, t0, t0, t0),w+1)/2 - 99*diff(F, t0, t0)._mul_trunc_(diff(F, t0, t0),w+1)/8 - 99*diff(F, t0)._mul_trunc_(diff(F, t0),w+1)._mul_trunc_(diff(F, t0, t0),w+1)/4 - 33*diff(F, t0).power_trunc(4,w+1)/8 - 891*diff(F, t0, t1).truncate(w+1)/2 - 891*diff(F, t1)._mul_trunc_(diff(F, t0),w+1)/2 + 1155*diff(F, t0).truncate(w+1)/32
    return Fs6(z=1)

def get_Fs24(F,w):
    Fs24 = 385*diff(F,t0)._mul_trunc_(diff(F,t0),w+1)/8 + 385*diff(F,t0,t0).truncate(w+1)/8 - 21*diff(F,t0)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t0,t0),w+1)/4 - 7*diff(F,t0)._mul_trunc_(diff(F,t0,t0),w+1)._mul_trunc_(diff(F,t0,t0),w+1) - 7*diff(F,t0,t0,t0,t0,t0).truncate(w+1)/12 - 147*diff(F,t0)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t1),w+1) - 637*diff(F,t0,t0,t1).truncate(w+1)/4 + 102*diff(F,t0)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t2),w+1) + 18*diff(F,t0,t1)._mul_trunc_(diff(F,t0,t1),w+1) - 7*diff(F,t0)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t0),w+1)/2 + 60*diff(F,t2)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t0),w+1) - 35*diff(F,t0)._mul_trunc_(diff(F,t0,t0,t0,t0),w+1)/12 - 21*diff(F,t0,t0)._mul_trunc_(diff(F,t0,t0,t0),w+1)/4 - 637*diff(F,t0)._mul_trunc_(diff(F,t0,t1),w+1)/2 + 36*diff(F,t1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t1),w+1) + 132*diff(F,t0)._mul_trunc_(diff(F,t0,t0,t2),w+1) + 30*diff(F,t2)._mul_trunc_(diff(F,t0,t0,t0),w+1) + 102*diff(F,t0,t2)._mul_trunc_(diff(F,t0,t0),w+1) + 18*diff(F,t0)._mul_trunc_(diff(F,t0,t1,t1),w+1) + 18*diff(F,t1)._mul_trunc_(diff(F,t0,t0,t1),w+1) - 720*diff(F,t0,t2)._mul_trunc_(diff(F,t2),w+1) - 1680*diff(F,t0)._mul_trunc_(diff(F,t0,t4),w+1) - 720*diff(F,t2,t2)._mul_trunc_(diff(F,t0),w+1) - 432*diff(F,t1,t2)._mul_trunc_(diff(F,t1),w+1) - 720*diff(F,t0,t2,t2).truncate(w+1) - 840*diff(F,t0,t0,t4).truncate(w+1) + 44*diff(F,t2,t0,t0,t0).truncate(w+1) + 20160*diff(F,t2,t4).truncate(w+1) - 216*diff(F,t1,t1,t2).truncate(w+1) - 147*diff(F,t1)._mul_trunc_(diff(F,t0,t0),w+1) + 9*diff(F,t1,t1,t0,t0).truncate(w+1) + 2520*diff(F,t2)._mul_trunc_(diff(F,t1),w+1) + 6720*diff(F,t0)._mul_trunc_(diff(F,t3),w+1) + 6720*diff(F,t0,t3).truncate(w+1) - 2835*diff(F,t2).truncate(w+1) + 2814*diff(F,t1,t2).truncate(w+1) - 332640*diff(F,t5).truncate(w+1)
    return Fs24(z=1)

def get_Fs33(F,w):
    Fs33 = ZZ(1)/2*(14400*diff(F,t3,t3).truncate(w+1) - 332640*diff(F,t5).truncate(w+1) - 1440*diff(F,t0,t1,t3).truncate(w+1) - 1440*diff(F,t1,t3)._mul_trunc_(diff(F,t0),w+1) - 1440*diff(F,t0,t3)._mul_trunc_(diff(F,t1),w+1) + 7020*diff(F,t0,t3).truncate(w+1) + 6720*diff(F,t0)._mul_trunc_(diff(F,t3),w+1) + 2160*diff(F,t1,t2).truncate(w+1) + 2160*diff(F,t2)._mul_trunc_(diff(F,t1),w+1) + 36*diff(F,t1,t1,t0,t0).truncate(w+1) + 36*diff(F,t1,t0)._mul_trunc_(diff(F,t1,t0),w+1) + 36*diff(F,t1,t1)._mul_trunc_(diff(F,t0,t0),w+1) + 72*diff(F,t1)._mul_trunc_(diff(F,t0,t0,t1),w+1) + 72*diff(F,t0)._mul_trunc_(diff(F,t0,t1,t1),w+1) + 36*diff(F,t1)._mul_trunc_(diff(F,t1),w+1)._mul_trunc_(diff(F,t0,t0),w+1) + 72*diff(F,t1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t1),w+1) + 36*diff(F,t1,t1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0),w+1) - 2400*diff(F,t2).truncate(w+1) - 165*diff(F,t0,t0,t1).truncate(w+1) - 165*diff(F,t1)._mul_trunc_(diff(F,t0,t0),w+1) - 150*diff(F,t1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0),w+1) - 315*diff(F,t0)._mul_trunc_(diff(F,t0,t1),w+1) + 725*diff(F,t0,t0).truncate(w+1)/16 + 175*diff(F,t0)._mul_trunc_(diff(F,t0),w+1)/4)
    return Fs33(z=1)

def get_Fs7(F,w):
    Fs7 = 17297280*diff(F,t7).truncate(w+1) - 50400*diff(F,t2,t3).truncate(w+1) - 50400*diff(F,t3)._mul_trunc_(diff(F,t2),w+1) - 90720*diff(F,t1,t4).truncate(w+1) - 90720*diff(F,t4)._mul_trunc_(diff(F,t1),w+1) - 332640*diff(F,t0,t5).truncate(w+1) - 332640*diff(F,t5)._mul_trunc_(diff(F,t0),w+1) + 468*diff(F,t1,t1,t1).truncate(w+1) + 1404*diff(F,t1)._mul_trunc_(diff(F,t1,t1),w+1) + 468*diff(F,t1).power_trunc(3,w+1) + 4680*diff(F,t0,t1,t2).truncate(w+1) + 4680*diff(F,t1)._mul_trunc_(diff(F,t0,t2),w+1) + 4680*diff(F,t2)._mul_trunc_(diff(F,t0,t1),w+1) + 4680*diff(F,t1,t2)._mul_trunc_(diff(F,t0),w+1) + 4680*diff(F,t2)._mul_trunc_(diff(F,t1),w+1)._mul_trunc_(diff(F,t0),w+1) + 5460*diff(F,t0,t0,t3).truncate(w+1) + 5460*diff(F,t3)._mul_trunc_(diff(F,t0,t0),w+1) + 10920*diff(F,t0)._mul_trunc_(diff(F,t0,t3),w+1) + 5460*diff(F,t3)._mul_trunc_(diff(F,t0),w+1).power_trunc(2,w+1) + 507780*diff(F,t4).truncate(w+1) - 10725*diff(F,t0,t2).truncate(w+1) - 10725*diff(F,t2)._mul_trunc_(diff(F,t0),w+1) - 5577*diff(F,t1,t1).truncate(w+1)/2 - 5577*diff(F,t1).power_trunc(2,w+1)/2 - 143*diff(F,t0,t0,t0,t1).truncate(w+1) - 143*diff(F,t1)._mul_trunc_(diff(F,t0,t0,t0),w+1) - 429*diff(F,t0,t1)._mul_trunc_(diff(F,t0,t0),w+1) - 429*diff(F,t0)._mul_trunc_(diff(F,t0,t0,t1),w+1) - 429*diff(F,t0).power_trunc(2,w+1)._mul_trunc_(diff(F,t0,t1),w+1) - 429*diff(F,t1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t0),w+1) - 143*diff(F,t1)._mul_trunc_(diff(F,t0),w+1).power_trunc(3,w+1) + 1001*diff(F,t0,t0,t0).truncate(w+1)/8 + 3003*diff(F,t0)._mul_trunc_(diff(F,t0,t0),w+1)/8 + 1001*diff(F,t0).power_trunc(3,w+1)/8 + 27027*diff(F,t1).truncate(w+1)/16 - ZZ(5005)/384
    return Fs7(z=1)

def get_Fs222(F,w):
    Fs222 = ZZ(1)/6*(1728*diff(F,t2,t2,t2).truncate(w+1) - 30240*diff(F,t2,t3).truncate(w+1) - 216*diff(F,t0,t0,t2,t2).truncate(w+1) - 432*diff(F,t0,t2).power_trunc(2,w+1) - 432*diff(F,t0)._mul_trunc_(diff(F,t0,t2,t2),w+1) + 864*diff(F,t0,t1,t2).truncate(w+1) + 864*diff(F,t1)._mul_trunc_(diff(F,t0,t2),w+1) + 864*diff(F,t1,t2)._mul_trunc_(diff(F,t0),w+1) + 1260*diff(F,t0,t0,t3).truncate(w+1) + 2520*diff(F,t0)._mul_trunc_(diff(F,t0,t3),w+1) + 9*diff(F,t0,t0,t0,t0,t2).truncate(w+1) + 36*diff(F,t0,t2)._mul_trunc_(diff(F,t0,t0,t0),w+1) + 36*diff(F,t0)._mul_trunc_(diff(F,t0,t0,t0,t2),w+1) + 36*diff(F,t0,t0)._mul_trunc_(diff(F,t0,t0,t2),w+1) + 72*diff(F,t0)._mul_trunc_(diff(F,t0,t2),w+1)._mul_trunc_(diff(F,t0,t0),w+1) + 36*diff(F,t0).power_trunc(2,w+1)._mul_trunc_(diff(F,t0,t0,t2),w+1) + 151200*diff(F,t4).truncate(w+1) - 576*diff(F,t1,t1).truncate(w+1) - 576*diff(F,t1).power_trunc(2,w+1) - 2628*diff(F,t0,t2).truncate(w+1) - 2520*diff(F,t2)._mul_trunc_(diff(F,t0),w+1) - 36*diff(F,t0,t0,t0,t1).truncate(w+1) - 108*diff(F,t0)._mul_trunc_(diff(F,t0,t0,t1),w+1) - 36*diff(F,t1)._mul_trunc_(diff(F,t0,t0,t0),w+1) - 72*diff(F,t0,t1)._mul_trunc_(diff(F,t0,t0),w+1) - 72*diff(F,t1)._mul_trunc_(diff(F,t0),w+1)._mul_trunc_(diff(F,t0,t0),w+1) - 72*diff(F,t0).power_trunc(2,w+1)._mul_trunc_(diff(F,t0,t1),w+1) - diff(F,t0,t0,t0,t0,t0,t0).truncate(w+1)/8 - 3*diff(F,t0)._mul_trunc_(diff(F,t0,t0,t0,t0,t0),w+1)/4 - 3*diff(F,t0,t0)._mul_trunc_(diff(F,t0,t0,t0,t0),w+1)/2 - 5*diff(F,t0,t0,t0).power_trunc(2,w+1)/4 - 6*diff(F,t0)._mul_trunc_(diff(F,t0,t0),w+1)._mul_trunc_(diff(F,t0,t0,t0),w+1) - diff(F,t0,t0).power_trunc(3,w+1) - 3*diff(F,t0).power_trunc(2,w+1)._mul_trunc_(diff(F,t0,t0,t0,t0),w+1)/2 - 3*diff(F,t0).power_trunc(2,w+1)._mul_trunc_(diff(F,t0,t0).power_trunc(2,w+1),w+1) - diff(F,t0).power_trunc(3,w+1)._mul_trunc_(diff(F,t0,t0,t0),w+1) + 63*diff(F,t0,t0,t0).truncate(w+1)/2 + 90*diff(F,t0)._mul_trunc_(diff(F,t0,t0),w+1) + 27*diff(F,t0).power_trunc(3,w+1) + 378*diff(F,t1).truncate(w+1) - ZZ(63)/20)
    return Fs222(z=1)

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
        F = self.F_series.get((),T.zero()) 
        if w > F_max_weight:
            tic = time.time()
#             time_est = time_for_F(w) - time_for_F(F_max_weight)
#             if time_est > 120:
#                 command = input(f"    Partition function update might take more than {float2time(time_est,2)}. Do you want to \
#                 continue? Print 'n' to abort, to continue press Enter.")
#                 if command == 'n':
#                     sys.exit("The computation was aborted by user due to time constraints.")
            if self.verbose: print(f"    Updating the partition function F from max_weight {F_max_weight} to {w}...")
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
        Fs = self.F_series.get(s_part,T.zero()) 
        if w > Fs_max_weight:
            tic = time.time()
#             time_est = self.times[s_part](w)
#             if time_est > 120:
#                 command = input(f"Fs function update for s = {s_part}, w = {w} might take more than {float2time(time_est,2)}.\
#                 Do you want to continue? Print 'n' to abort, to continue press Enter..")              
#                 if command == 'n':
#                     sys.exit("The computation was aborted by user due to time constraints.")
            if self.verbose: print(f"Updating Fs function for s = {s_part} from max_weight {Fs_max_weight} to {w}...")
            F = self.partition_function(w+self.shifts[s_part])
            Fs = self.AC_formulae[s_part](F,w)
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
