from sage.all import *

R = PolynomialRing(QQ, ['t%d' % i for i in range(31)])
S = PolynomialRing(QQ, ['b%d' % i for i in range(1,31)])

def float2time(t,precision=0):
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
    return ZZ(n).multifactorial(2)

def approx(x):
    return "{0:.4f}".format(float(N(x)))

# c_f(d) returns the renormalization coefficient 2^(d+1)/(d-1)!
def c_d(stratum):
    return sum(stratum)/2+len(stratum)

def c_f(dim_or_stratum):
    if type(dim_or_stratum) == list:
        d = c_d(dim_or_stratum)
    return 2**(d+1)/factorial(d-1)