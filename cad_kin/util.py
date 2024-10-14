import numpy as np



def decompose(a,b):
    '''
    a: 2D vector
    b: 2D vector

    returns decomposes b into components orthogonal and inline with a in two dimensions
    '''
    a_unit = a/np.linalg.norm(a)
    a_unit = a_unit[:,None]
    b = b[:,None]
    proj = np.matmul(np.matmul(a_unit,np.transpose(a_unit)),b)
    orth = b-proj
    return orth,proj