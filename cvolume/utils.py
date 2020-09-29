from sage.all import PolynomialRing, ZZ, QQ, N, exp, factorial

R = PolynomialRing(QQ, ['t%d' % i for i in range(31)])
S = PolynomialRing(QQ, ['b%d' % i for i in range(1,31)])

def float2time(t,precision=0):
    '''
    Takes float and precision (int) and return a time string in minutes and seconds.
    '''
    minutes, seconds = divmod(int(t*10**precision), 60*10**precision)
    if minutes: return '%s min %.*f s' % (int(minutes), int(precision), int(seconds)/10**precision)
    else: return '%.*f s' % (int(precision),int(seconds)/10**precision)

time_for_F = lambda w : N(exp(0.39*w-6.7))
time_for_Fs2 = lambda w : N(exp(0.38*w-5.2))
time_for_Fs3 = lambda w : N(exp(0.38*w-5))
time_for_Fs4 = lambda w : N(exp(0.71*w-6.4))
time_for_Fs22 = lambda w : N(exp(0.63*w-4.8))
time_for_Fs5 = lambda w : N(exp(0.68*w-5.5))
time_for_Fs23 = lambda w : N(exp(0.75*w-4.6))
time_for_Fs6 = lambda w : N(exp(0.885*w-3.8))
time_for_Fs24 = lambda w : N(exp(0.96*w-3.6))
time_for_Fs33 = lambda w : N(exp(0.745*w-3.8))
time_for_Fs7 = lambda w : N(exp(0.96*w-3.6))
time_for_Fs222 = lambda w : N(exp(w-2.25))

def dfactorial(n):
    '''
    Return double-factorial of integer.
    '''
    return ZZ(n).multifactorial(2)

def c_d(stratum):
    '''
    Return the dimension of a stratum.
    '''
    return sum(stratum)/ZZ(2)+len(stratum)

def c_f(dim_or_stratum):
    '''
    Return the renormalization coefficient for stratum: 2^(d+1)/(d-1)! where d is dimension. Input can be a dimension or a stratum itself.
    '''
    if type(dim_or_stratum) == list:
        d = c_d(dim_or_stratum)
    else:
        d = dim_or_stratum
    return ZZ(2**(d+1))/factorial(d-1)

def coeff_2(g,n,M): 
    '''
    Return power-of-2 renormalizing coefficient for local polynomials.
    '''
    return ZZ(1)/2**(5*g-6+2*n-2*M)

