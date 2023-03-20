#pragma once


#include <stdio.h>
#include <stdlib.h>

typedef struct convolutional_code_s  convolutional_code;

/* output_bits is the number of bits the encoder output, 1/r 
 * poly array specifies the code poly nominal. maximum constrained 
 * length is 32
 * 
 * Suppose the generator is g(D) = a_{n}D^n + a_{n-1}D^{n-1] + a_0 
 * the polynominal is then |a_{n}|a_{n-1}|.... |a_0|
 * 
 * for systematic codes, the first one is the denominator, so a_0 should be 1
 * for poly[0]
 * 
 */

typedef struct{
    unsigned  nb_out_bits;
    unsigned *poly;
}convolutional_code_coefs;
    
convolutional_code* convoluational_code_construct( convoluational_code_coefs *p_coefs);


void convolutional_code_destroy(convolutional_code *p_code);

/* input is the octets of input bits, so is the output, 
 * the sz specitifies the number of bytes of the input array 
 * and the pre-allocated bytes for the output,  the nb_bits specifies how
 * many bits need to encodied
 * 
 * the rsc version is the systematic encoder, information bits will be the 
 * first output, followed by the parity bits. 
 *
 */
unsigned encoding_bits(convolutional_code *p_code
                      ,unsigned char      *input
                      ,size_t              input_sz
                      ,unsigned char      *output
                      ,size_t              output_sz
                      ,unsigned            nb_bits
                      );

unsigned encoding_bits_rsc(convolutional_code *p_code
                          ,unsigned char      *input
                          ,size_t              input_sz
                          ,unsigned char      *output
                          ,size_t              output_sz
                          ,unsigned            nb_bits
                          );





   
