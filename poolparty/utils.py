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