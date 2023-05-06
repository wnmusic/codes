import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np
import pandas
from pygnuplot import gnuplot

def get_trellis(taps, b_rsc):
    
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
    trellis = []
    for i in range(1<<order):
        s0,s1,ibit,obits0, obits1 = get_trellis_node(cc, i)
        print("s: {:02b}, s0: {:02b}, s1: {:02b}, i: {:01b}, o0: {:02b}, o1: {:02b}".format(i, s0, s1, ibit, obits0, obits1));
        trellis.append((i, s0, s1, ibit, obits0, obits1))

    convolutional_code_destroy(cc)
    return trellis


if __name__ == '__main__':
    if (len(sys.argv) < 3):
        print("Usage: plot_trellis  '[tap0,tap1]', b_rsc out.png")

    taps = eval(sys.argv[1])
    b_rsc = int(sys.argv[2])
    if (len(sys.argv) > 3):
        outfile = sys.argv[3]

    trans = get_trellis(taps, b_rsc)
    num_state = len(trans)
    
    gp = gnuplot.Gnuplot()
    gp.set('term png enhanced font "arial,15"')
    gp.set('output "{}"'.format(outfile))
    gp.unset('tics')
    gp.set('xrange [-0.2:1.2] writeback')
    gp.set('yrange [-0.2:{}] writeback'.format(num_state-0.8))


    states = np.asarray(2 * [[i for i in range(num_state)]])
    df = pandas.DataFrame(states)
    cmd = ['u 1:{} w p pt 6 ps 4 notitle'.format(i) for i in range(2,num_state+2)]

    set_cmd = []
    for t in trans:
        s = t[0]
        prev_s0 = t[1]
        prev_s1 = t[2]
        ibit = t[3]
        o0 = t[4]
        o1 = t[5]
        set_cmd += ['label "{}" at first 0,{} center front'.format(s,s)]
        set_cmd += ['label "{}" at first 1,{} center front'.format(s,s)]
        set_cmd += ['arrow from first 1,{} to first 0,{}  backhead front'.format(s, prev_s0)]
        if s == num_state - 1:
            pos = s - 0.2
        else:
            pos = s
        
        set_cmd += ['label "{:01b}/{:02b}" at first 0.8,{} '.format(ibit, o0, pos+0.2*( 1 if prev_s0 > prev_s1 else 0))]
        set_cmd += ['arrow from first 1,{} to first 0,{}  backhead front'.format(s, prev_s1)]
        set_cmd += ['label "{:01b}/{:02b}" at first 0.8,{} '.format(ibit, o1, pos+0.2*( 1 if prev_s1 > prev_s0 else 0))]        

    gp.set(*set_cmd)
    gp.plot_data(df, *cmd)

        
