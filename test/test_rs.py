import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np



def test_rs_encoding(n, k):
    rs = rs_code_construct(n, k)
    info = np.zeros(k, dtype=np.uint8);
    info[1] = 2;
    info[0] = 7;


    d, code = rs_encode(rs, info, n)
    print(d, code)

    d, rec = rs_decode_ge(rs, code, k)
    print(d, rec)

    
    #d, genpoly = rs_get_genpoly(rs, n-k+1);
    #print(d, genpoly[0:d])
    
    #d, code = rs_encode_sys(rs, info, n)
    #print(d, code)
    
    rs_code_destroy(rs)
    

if __name__ == '__main__':
    input(os.getpid())
    test_rs_encoding(7, 3)
