import sys
import os
sys.path.append('../build/test/')
from pycode import *
import pytest
import math
import numpy as np



@pytest.mark.parametrize("taps, b_rsc", [([5, 7], 0)
                                         ,([5, 7], 1)])
def test_encoding(taps, b_rsc):
    poly = uint32Array(len(taps))
    rate = len(taps)
    for i,t in enumerate(taps):
        poly[i] = t
        
    coefs = convolutional_code_coefs()
    coefs.nb_out_bits = len(taps)
    coefs.poly = poly
    coefs.b_rsc = b_rsc
    coefs.dec_delay_factor=5
    
    cc = convolutional_code_construct(coefs)

    n_in_bytes = 1
    n_in_bits = 8
    in_bytes = np.zeros(n_in_bytes, dtype=np.ubyte)
    out_sz = rate*n_in_bytes + 1

    n_out_bits,out_bytes = encode(cc, in_bytes, out_sz, n_in_bits)
    print('test all zeros input, #of output bits {}'.format(n_out_bits))
    for i in range((n_out_bits + 7)//8):
        print("{:08b}".format(out_bytes[i]))                


    n_in_bits = 1
    in_bytes[0] = (1<<7)
    print('test imp response input {}'.format(in_bytes))

    n_out_bits,out_bytes = encode(cc, in_bytes, out_sz, n_in_bits)
    print('# of output bits {}'.format(n_out_bits))        
    for i in range((n_out_bits + 7)//8):
        print("{0:08b}".format(out_bytes[i]))        
        
    convolutional_code_destroy(cc)


def print_trellis(taps, b_rsc):
    
    poly = uint32Array(len(taps))
    rate = len(taps)
    for i,t in enumerate(taps):
        poly[i] = t

    order = (int)(np.floor(np.log2(max(taps))))
    coefs = convolutional_code_coefs()
    coefs.nb_out_bits = len(taps)
    coefs.poly = poly
    coefs.b_rsc = b_rsc
    coefs.dec_delay_factor=5
    
    cc = convolutional_code_construct(coefs)

    for i in range(1<<order):
        s0,s1,ibit,obits0, obits1 = get_trellis_node(cc, i)
        print("s: {:02b}, s0: {:02b}, s1: {:02b}, i: {:01b}, o0: {:02b}, o1: {:02b}".format(i, s0, s1, ibit, obits0, obits1));

    convolutional_code_destroy(cc)
    
if __name__ == '__main__':
    input(os.getpid())
    test_encoding([5, 7], 1)
    print_trellis([5, 7], 1)
