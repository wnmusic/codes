import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np



def test_rs255_encoding(k):
    rs = rs255_encoder_construct(k)
    info = np.zeros(k, dtype=np.uint8);
    info[2] = 1;

    n, genpoly = rs255_encoder_genpoly(rs, 40);
    print(n, genpoly)

    #n, code = rs255_encode(rs, info, 255)
    #print(n, code)

    rs255_encoder_destroy(rs)
    


if __name__ == '__main__':
    input(os.getpid())
    test_rs255_encoding(233)
