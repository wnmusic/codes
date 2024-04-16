import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np



def test_rs_encoding(n, k):
    rs = rs_code_construct(n, k)
    info = np.random.randint(0, n, k, dtype=np.uint8);
    
    max_t = (n-k+1)//2

    err = np.zeros(n, dtype=np.uint8);

    for i in range(max_t):
        pos = np.random.randint(0, n)
        err[pos] = np.random.randint(0, n)
        
    _, code = rs_encode(rs, info, n)

    code = np.bitwise_xor(code, err)

    _, rec = rs_decode_ge(rs, code, k)

    assert np.any(info - rec)==False, "something defintely wrong\n info = {} \n rec = {}".format(info, rec)

    
    #d, genpoly = rs_get_genpoly(rs, n-k+1);
    #print(d, genpoly[0:d])
    
    #d, code = rs_encode_sys(rs, info, n)
    #print(d, code)
    
    rs_code_destroy(rs)
    

if __name__ == '__main__':
    input(os.getpid())
    test_rs_encoding(255, 233)
