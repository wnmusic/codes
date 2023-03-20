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

%apply (unsigned char *IN_ARRAY1, int DIM1) {(unsigned char *input, int input_sz)};
%apply (unsigned char *ARGOUT_ARRAY1, int DIM1) {(unsigned char *output, int output_sz)};

%include "conv.h"

