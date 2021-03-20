#Copyright 2021 Lixian WANG. All Rights Reserved.
import warnings
import functools
from itertools import chain as ic
import numpy as np

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)  # turn off filter
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning,
                      stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)
    return new_func

def flattenList(nlist_input):
    nlist = nlist_input
    if not isinstance(nlist,list):
        raise ValueError(f'this is not a list')
    try:
        while True:
            nlist = list(ic(*nlist))
    except TypeError:
        pass
    return nlist

def div(xlist,number):
    if not isinstance(xlist,(list,np.ndarray,np.generic)):
        raise ValueError('xlist is not a list or np.ndarray object')
    if not isinstance(number,(int,float)):
        raise ValueError('number is not an int or a float object')
    return list(map(lambda x: x/number,xlist))
