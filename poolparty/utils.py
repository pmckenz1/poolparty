from numba import njit
import numpy as np
from tqdm import tqdm_notebook

def get_read(arr,
             ncutsites,
            ):
    '''
    arr is shape aliquots x gamete indexes
    this fuction returns the aliquot number, the gamete index, and the cut site index (from 0 to `ncutsites`)
    '''
    alinum = np.random.randint(arr.shape[0])
    colnum = np.random.randint(arr.shape[1])
    gamidx = arr[alinum,colnum]
    cutsiteidx = np.random.randint(ncutsites)
    return(alinum, gamidx, cutsiteidx)


@njit
def convert_haplo(n,s):
    s_ = s
    for i in range(n):
        s_ = 1-s_
    return(s_)


# because scipy.stats.binom is slow...
def poisson_pmf(k, lam):
    try:
        x = (lam**k) * np.exp(-1*lam) / np.math.factorial(k)
        return(x)
    except:
        return( np.nan )

@njit
def binom_jit(n,k,p,com):
    return(com*p**k*(1-p)**(n-k))


# full binomial function!
def binomial(n,k,p):
    com=binom(n,k)
    return(binom_jit(n,k,p,com))    


# get the combinations of different crossovers from n_1 and ali_size that would yield a given n_2
# (that is, the different ways hap1s can switch to hap0 and hap0s can switch to hap1 to give us n_2)
@njit
def get_combs(n_1,ali_size,n_2):
    zs = ali_size - n_1
    arr_ = np.zeros((n_1+1,zs+1),dtype=np.int64)
    for row in range(arr_.shape[0]):
        for col in range(arr_.shape[1]):
            arr_[row,col] = n_1 - row + col
    # top row is the number of n_1 that change to zero, bottom is number of zero that change to 1
    return(np.vstack(np.where(arr_ == n_2)).T)


# calculate the probability of switching to an n_2 value given an n_1 value, p_1, p_2, and aliquot size
def calc_n_2_prob(n_1, n_2, p_1, p_2, ali_size):
    sumprob_n_2 = 0.0
    for i in get_combs(n_1, ali_size, n_2):
        nstay = n_1 - i[0] # nstay is the number of n_1 that stay n_1... i[0] is the number that change!
        nchange = i[1]
        # prob that nstay ones _do_ stay at one
        p_nstay = binomial(k=nstay, n=n_1, p=(1-p_2)*(1-p_1)+(1-p_2)*p_1*(1-p_2)+p_2*p_1*p_2)
        # prob that nchange zeros _do_ change to one
        p_nchange = binomial(k=nchange, n=ali_size-n_1, p=p_2*(1-p_1)+p_2*p_1*(1-p_2)+(1-p_2)*p_1*p_2)

        sumprob_n_2 += p_nstay * p_nchange
    return(sumprob_n_2)