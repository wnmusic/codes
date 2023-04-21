%module pycode
%{
#define SWIG_FILE_WITH_INIT
#include "conv.h"
%}

%include "numpy.i"
%include "carrays.i"
%array_class(unsigned, uint32Array)
%init %{
import_array();
%}

%apply (unsigned char *IN_ARRAY1, int DIM1) {(unsigned char *in, int in_sz)};
%apply (unsigned char *ARGOUT_ARRAY1, int DIM1) {(unsigned char *out, int out_sz)};

%include "typemaps.i"
%apply (unsigned *OUTPUT) { unsigned *p0, unsigned *p1, unsigned  *p_in_bit, unsigned  *p_out0_bits, unsigned *p_out1_bits};

%include "conv.h"

