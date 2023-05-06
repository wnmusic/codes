import sys
import os
sys.path.append('../build/test/')
from pycode import *
import pytest
import math
import numpy as np



@pytest.mark.parametrize("taps, b_rsc", [([5, 7], 0)
                                         ,([5, 7], 1)])
def test_viterbi(taps, b_rsc, ebn0_dB):
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

    n_in_bytes = 2
    n_in_bits = 16
    in_bytes = np.zeros(n_in_bytes, dtype=np.ubyte)
    out_sz = rate*n_in_bytes + 1

    in_bytes[0] = 255
    in_bytes[1] = 33
    

    print('information bits {}'.format(n_in_bits))
    for i in range((n_in_bits + 7)//8):
        print("{:08b}".format(in_bytes[i]))
    
    n_out_bits,out_bytes = encode(cc, in_bytes, out_sz, n_in_bits)
    print('encoded bits {}'.format(n_out_bits))
    for i in range((n_out_bits + 7)//8):
        print("{:08b}".format(out_bytes[i]))

    N0_2 = rate * np.power(10.0, -ebn0_dB * 0.1) * 0.5
    np.random.seed(10)

    xmit_symbols = []
    for o in out_bytes:
        for i in range(8):
            xmit_symbols +=[ 1.0 - ((o>>i) & 1)*2.0 ]

    print('transmit: {}'.format(xmit_symbols[:n_out_bits]))
    recv_symbols = np.array(xmit_symbols, dtype=np.float32) + np.float32(N0_2 * np.random.randn(len(xmit_symbols)))
    print('recv: {}'.format(recv_symbols[:n_out_bits]))
        
    n_dec_bits, dec_bytes = viterbi_decode(cc, recv_symbols, out_sz, n_out_bits)
    print('decoded bits {}:'.format(n_dec_bits))
    for i in range((n_dec_bits + 7)//8):
        print("{:08b}".format(dec_bytes[i]))
        
    convolutional_code_destroy(cc)

if __name__ == '__main__':
    input(os.getpid())
    test_viterbi([7, 5], 1, 6.0)

