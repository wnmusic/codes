import os, sys
import numpy as np
import math
import struct
# the trace head is 64 of struct
#typedef struct{
#   char name[64];
#    unsigned ele_sz; /* sizeof each element */
#    unsigned nb_ele; /* how many elements per frame */
#    char ctype[16]; /* float or int ... */
#}trace_info;
import ctypes

def read_trace(fname):
    with open(fname, mode='rb') as f:
        buf = f.read()

    offset = struct.calcsize('64sII16s')*64
    header_buf = buf[:offset]
    hiter = struct.iter_unpack('64sII16s', header_buf)
    header = []
    for i in hiter:
        name = ctypes.create_string_buffer(i[0]).value
        if name == b'':
            break;
        ele_sz = i[1]
        ele_nb = i[2]
        ctype = ctypes.create_string_buffer(i[3]).value
        if ctype == b'float' and ele_sz == 4:
            dtype = np.float32
        elif (ctype == b'int' or ctype == b'int32_t') and ele_sz == 4:
            dtype = np.int32
        elif (ctype == b'short' or ctype == b'int16_t') and ele_sz == 2:
            dtype = np.int16
        elif (ctype == b'unsigned short' or ctype == b'uint16_t')  and ele_sz == 2:
            dtype = np.uint16
        header += [(name, dtype, ele_sz, ele_nb)]

    ret = {}
    while(offset < len(buf)):
        for h in header:
            tag = h[0]
            dt = h[1]
            dim = h[3]
            a = np.frombuffer(buf, dtype=dt, count=dim, offset = offset)
            if tag not in ret:
                ret[tag] = a
            else:
                ret[tag] = np.vstack((ret[tag], a))
            offset = offset + dim * h[2]
    return ret

if __name__ == '__main__':
    fname = sys.argv[1]
    print(read_trace(fname))

    
    

