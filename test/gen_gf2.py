import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np

def generate_primitive_mapping(deg):
    if (deg > 8):
        print("Not supported yet")
        return
    
    n, prime = find_prime_poly(deg, 256)
    a = "{} prime polynomials of degree {} are: \n".format(n, deg)
    for i in range(n):
        a += "{}) : {}\n".format(i, bin(prime[i]))
    a += "choose one from above:\n"
    
    sel = input(a)
    p = int(prime[int(sel)]);
    
    la = [1]
    alpha = find_primitive_element(p);
    la += [alpha]
    alpha_1 = alpha
    for i in range(2,(1<<deg)-1):
        alpha_1 = gf2_8_mult_mod(alpha_1, alpha, p)
        la += [alpha_1]

    cstr = "static const uint8_t alphas_%d[] = { \n"%deg
    for (i, e) in enumerate(la):
        cstr += "%d, " % e
        if (i % 16 == 15):
            cstr += '\n'
    cstr += '\n};'

    print(cstr)

    # genreate the reverse index of the alphas
    cstr = "static const uint8_t index_of_alphas_%d[] = { 255, \n"%deg
    for i in range(1, 1<<deg):
        cstr += "%d, " % la.index(i)
        if (i % 16 == 15):
            cstr += '\n'
    cstr += '\n};'
    print(cstr)


if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print("usage: gen_gf2.py deg (at most 8)")
        sys.exit(1)
    
    generate_primitive_mapping(int(sys.argv[1]))
    sys.exit(0)


