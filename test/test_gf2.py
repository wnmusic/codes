import sys
import os
sys.path.append('../build/test/')
from pycode import *
import math
import numpy as np


def test_poly_mult():
    
    f = np.array([2, 34], dtype=np.uint8)
    g = np.array([88, 48, 99], dtype=np.uint8)
    n, fg = gf2m_poly_mult(f, g, 5)
    print(n, fg)

def test_gf2m_poly_div():
    
    f = np.array([29, 124], dtype=np.uint8)
    g = np.array([88, 48, 99, 246, 2, 1], dtype=np.uint8)
    q, dq, dr = gf2m_poly_div(g, f, 5)
    print(q, dq, dr, g)

def test_gf2m_ge():
    A = np.array([[88, 48, 99, 125],
                   [246, 2, 1, 5],
                   [137, 72, 9, 4],
                  [43, 47, 234, 6]
                  ], dtype=np.uint8)
    b = np.array([43, 89, 101, 27], dtype=np.uint8)
    gf2m_gaussian_elemination(A, b)
    print(A, b)
    

if __name__ == '__main__':
    input(os.getpid())
    #test_poly_mult()
    #test_gf2m_poly_div()
    test_gf2m_ge()

    
    #nd = 8;
    #n, prime = find_prime_poly(nd, 256)
    #print("{} prime polynomials".format(n))
    #print("using the first prime polynomial {}".format(bin(prime[0])))
    #
    #p = int(prime[0])
    #
    #la = [1]
    #alpha = find_primitive_element(p);
    #la += [alpha]
    #alpha_1 = alpha
    #for i in range(2,(1<<nd)-1):
    #    alpha_1 = gf2_8_mult_mod(alpha_1, alpha, p)
    #    la += [alpha_1]
    #
    #
    #print(la)
    #print(len(la))
    #print(sorted(la))


                
    
