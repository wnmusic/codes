import sys
import os
sys.path.append('../release/test/')
from pycode import *
import math
import numpy as np
from pygnuplot import gnuplot
import pandas

def sim_one_bpsk_frame(ebn0_dB, nbytes):
    
    in_bytes = np.random.randint(0, 256, nbytes, dtype=np.ubyte)
    # this is raw message
    in_bits = np.unpackbits(in_bytes, bitorder='little')
    xmit_symbols = 1.0 - 2.0 * in_bits
    N0_2 = np.power(10.0, -ebn0_dB * 0.1) * 0.5
    
    recv_symbols = np.array(xmit_symbols, dtype=np.float32) + np.float32(np.sqrt(N0_2) * np.random.randn(len(xmit_symbols)))
    dec_bits = recv_symbols < 0
    
    return  np.count_nonzero(dec_bits - in_bits)

        
def sim_one_coded_frame(cc, rate, ebn0_dB, nbytes):
    
    n_in_bits = nbytes * 8
    out_sz = rate*nbytes + 1

    in_bytes = np.random.randint(0, 256, nbytes, dtype=np.ubyte)
    n_out_bits,out_bytes = encode(cc, in_bytes, out_sz, n_in_bits)
    N0_2 = rate * np.power(10.0, -ebn0_dB * 0.1) * 0.5

    xmit_symbols = 1.0 - 2.0*np.unpackbits(out_bytes, bitorder='little')
    
    recv_symbols = np.array(xmit_symbols, dtype=np.float32) + np.float32(np.sqrt(N0_2) * np.random.randn(len(xmit_symbols)))
    n_dec_bits, dec_bytes = viterbi_decode(cc, recv_symbols, out_sz, n_out_bits)
    
    dec_bits = np.unpackbits(dec_bytes[:nbytes], bitorder='little')
    in_bits = np.unpackbits(in_bytes, bitorder='little')
    err_bits = dec_bits - in_bits
    
    return  np.count_nonzero(err_bits)



if __name__ == '__main__':
    #input(os.getpid())
    np.random.seed(10)
    frame_len = 128

    EbN0_list = [7, 8, 9, 10, 11]
    num_frames_list = [1000, 1000, 10000, 50000, 100000]
    result = []
    
    for EbN0, num_frames in zip(EbN0_list, num_frames_list):
        blk_err = 0
        bit_err = 0
        
        for i in range(num_frames):
            ret = sim_one_bpsk_frame(EbN0, frame_len)
            bit_err += ret
            blk_err += ret > 0
            
        err = (EbN0, bit_err/num_frames/frame_len/8,  blk_err/num_frames)
        print("Eb/N0 {}dB bit error rate: {}, blk error rate {}".format(err[0], err[1], err[2]))
        result += [err]

    df = pandas.DataFrame(result, columns=['Eb/No', 'bit error rate', 'block error rate'])
    df.to_csv('uncoded_res.csv', sep=' ', index_label='#')

    
    taps = [5,7]
    poly = uint32Array(len(taps))
    rate = len(taps)
    for i,t in enumerate(taps):
        poly[i] = t
        
    coefs = convolutional_code_coefs()
    coefs.nb_out_bits = len(taps)
    coefs.poly = poly
    coefs.b_rsc = 0
    coefs.dec_delay_factor=10
    
    cc = convolutional_code_construct(coefs)
    
    num_frames_list = [500, 5000, 50000, 100000]
    EbN0_list = [4, 5, 6, 7]
    result = []
    
    for EbN0, num_frames in zip(EbN0_list, num_frames_list):
        
        blk_err = 0
        bit_err = 0
        
        for i in range(num_frames):
            ret = sim_one_coded_frame(cc, rate, EbN0, frame_len)
            bit_err += ret
            blk_err += ret > 0
            
        err = (EbN0, bit_err/num_frames/frame_len/8,  blk_err/num_frames)
        print("Eb/N0 {}dB bit error rate: {}, blk error rate {}".format(err[0], err[1], err[2]))
        result += [err]

    convolutional_code_destroy(cc)
    df = pandas.DataFrame(result, columns=['Eb/No', 'bit error rate', 'block error rate'])
    df.to_csv('cc_res.csv', sep=' ', index_label='#')


