import sys
import os
sys.path.append('../build/test/')
from pycode import *

input(os.getpid())
v = 2
poly = uint32Array(v)
poly[0] = 7
poly[1] = 5
coefs = convolutional_code_coefs()
coefs.nb_out_bits = 2
coefs.poly = poly
cc = convolutional_code_construct(coefs)
convolutional_code_destroy(cc)


