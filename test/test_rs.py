import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np


def test_rs_encoding(n, k):
    rs = rs_code_construct(n, k, RS_CODE_POLYVAL)
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

def test_rs_BM(n, k):
    rs = rs_code_construct(n, k, RS_CODE_CYCLIC)
    #d, genpoly = rs_get_genpoly(rs, n-k+1);
    #print(d, genpoly[0:d])
        
    info = np.random.randint(0, n, k, dtype=np.uint8);
    #info = np.array([6, 2, 0], dtype=np.uint8)
    max_t = (n-k+1)//2

    err = np.zeros(n, dtype=np.uint8);    
    for i in range(max_t):
        pos = np.random.randint(0, n)
        err[pos] = np.random.randint(0, n)

    #err = np.array([0, 0, 0, 4, 0, 3, 0], dtype=np.uint8)        

    _, code = rs_encode(rs, info, n)
    #print(info, code)
    code = np.bitwise_xor(code, err)
    #print(err, code)

    _, rec = rs_decode_BM(rs, code, k)

    assert np.any(info - rec)==False, "something defintely wrong\n diff = {}".format(np.bitwise_xor(info, rec))



def test_rs_all_erased(n, k):
    rs = rs_code_construct(n, k, RS_CODE_CYCLIC)
    #d, genpoly = rs_get_genpoly(rs, n-k+1);
    #print(d, genpoly[0:d])
        
    info = np.random.randint(0, n, k, dtype=np.uint8);
    #info = np.array([1, 6, 0], dtype=np.uint8)

    era = []
    for i in range(n-k):
        pos = np.random.randint(0, n)
        if (pos not in era):
            era += [pos];
    era = np.array(era, dtype=np.uint8);
    #era = np.array([6,3,2], dtype=np.uint8);

    _, code = rs_encode(rs, info, n)
    print(info,code, era)
    _, rec = rs_decode_BM_with_erasure(rs, code, era, k)

    assert np.any(info - rec)==False, "something defintely wrong\n diff = {}".format(np.bitwise_xor(info, rec))


def test_rs_with_erasure(n, k, max_nb_era):
    rs = rs_code_construct(n, k, RS_CODE_CYCLIC)
    #d, genpoly = rs_get_genpoly(rs, n-k+1);
    #print(d, genpoly[0:d])
        
    info = np.random.randint(0, n, k, dtype=np.uint8);
    #info = np.array([1, 6, 0], dtype=np.uint8)

    era = []
    for i in range(max_nb_era):
        pos = np.random.randint(0, n)
        if (pos not in era):
            era += [pos];
    era = np.array(era, dtype=np.uint8);
    #era = np.array([6,3,2], dtype=np.uint8);


    
    nb_era = np.count_nonzero(era);
    max_t = (n-k - nb_era)//2
    err = np.zeros(n, dtype=np.uint8);
    for i in range(max_t):
        pos = np.random.randint(0, n)
        if (pos not in era):
            err[pos] = np.random.randint(0, n)
    
    _, code = rs_encode(rs, info, n)
    print(info,code,era,err)
    _, rec = rs_decode_BM_with_erasure(rs, code, era, k)

    assert np.any(info - rec)==False, "something defintely wrong\n diff = {}".format(np.bitwise_xor(info, rec))
    
    
if __name__ == '__main__':
    input(os.getpid())
    #test_rs_encoding(7, 3)
    #test_rs_BM(7, 3)
    #test_rs_BM(255, 223)
    #test_rs_all_erased(7,3)
    #test_rs_with_erasure(7,3,2)
    #test_rs_with_erasure(15,7,1)
    test_rs_with_erasure(63,13,10)
    test_rs_with_erasure(255,223,5)
