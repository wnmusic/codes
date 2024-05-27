%module pycode
%{
#define SWIG_FILE_WITH_INIT
#include "conv.h"
#include "gf2.h"
#include "rs.h"
%}

%include "stdint.i"
%include "numpy.i"
%include "carrays.i"


%array_class(unsigned, uint32Array)
%init %{
import_array();
%}


%apply (unsigned char *IN_ARRAY1, int DIM1) {(unsigned char *in, int in_sz)};
%apply (unsigned char *IN_ARRAY1, int DIM1) {(uint8_t *f, int df)};
%apply (unsigned char *IN_ARRAY1, int DIM1) {(uint8_t *g, int dg)};
%apply (unsigned char *IN_ARRAY1, int DIM1) {(unsigned char *erasures, int nb_era)};

%apply (unsigned char *INPLACE_ARRAY1, int DIM1) {(uint8_t *fr, int dfr)};
%apply (unsigned char *ARGOUT_ARRAY1, int DIM1) {(unsigned char *out, int out_sz)};
%apply (unsigned char *ARGOUT_ARRAY1, int DIM1) {(uint8_t *out, int out_sz)};
%apply (unsigned *ARGOUT_ARRAY1, int DIM1) {(unsigned *out, int out_sz)};

%apply (float *IN_ARRAY1, int DIM1) {(float *input, int input_sz)};

%apply (unsigned char *INPLACE_ARRAY1, int DIM1) {(uint8_t *b, int dim)};

%apply (unsigned char *INPLACE_ARRAY2, int DIM1, int DIM2) {(uint8_t *A, int dim_A1, int dim_A2)};

%include "typemaps.i"
%apply (unsigned *OUTPUT) { unsigned *p0, unsigned *p1, unsigned  *p_in_bit, unsigned  *p_out0_bits, unsigned *p_out1_bits};
%apply (int *OUTPUT) { int *dq, int *dr};

%include "conv.h"
%include "gf2.h"
%include "rs.h"
