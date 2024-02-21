import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np



def test_rs_encoding(n, k):
    rs = rs_encoder_construct(n, k)
    info = np.zeros(k, dtype=np.uint8);
    info[2] = 1;

    n, genpoly = rs_get_genpoly(rs, n-k+1);
    print(n, genpoly)

    #n, code = rs_encode(rs, info, 255)
    #print(n, code)

    rs_encoder_destroy(rs)
    


if __name__ == '__main__':
    input(os.getpid())
    test_rs_encoding(7, 3)
