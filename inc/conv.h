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
    /* this specifies the rate of the code */
    unsigned  nb_out_bits;
    /* this is the taps */
    unsigned *poly;
    /* recursive systemaic code */
    int b_rsc;
    /* this times K is the decoder delay, -1 means maximum delay */
    unsigned  dec_delay_factor; 
}convolutional_code_coefs;
    
convolutional_code* convolutional_code_construct( convolutional_code_coefs *p_coefs);


void convolutional_code_destroy(convolutional_code *p_code);

void get_trellis_node(convolutional_code *p_code
                     ,unsigned node
                     ,unsigned *p0
                     ,unsigned *p1
                     ,unsigned      *p_in_bit
                     ,unsigned      *p_out0_bits
                     ,unsigned      *p_out1_bits
                     );

/* input is the octets of input bits, so is the output, 
 * the sz specitifies the number of bytes of the input array 
 * and the pre-allocated bytes for the output,  the nb_bits specifies how
 * many bits need to encodied
 * 
 * the rsc version is the systematic encoder, information bits will be the 
 * first output, followed by the parity bits. 

 * this encodes the input bits and terminate the block, so the number
 * of output bits should be (nb_bits + K)*v, where K is the constraint length 
 * v is the 1/rate, the output array should have enough space for the output
 *
 *
 */
unsigned encode(convolutional_code *p_code
               ,unsigned char      *in
               ,int              in_sz
               ,unsigned char      *out
               ,int              out_sz
               ,unsigned            nb_bits
               );






   
