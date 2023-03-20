#include "conv.h"
#include <string.h>



struct convolutional_code_s
{
    convolutional_code_coefs *p_coefs;
    unsigned *states;
    unsigned order;
};

convolutional_code* convolutional_code_construct( convolutional_code_coefs *p_coefs)
{
    convolutional_code *p_code = calloc(1, sizeof(convolutional_code));
    unsigned n_max = 0,  n;
    p_code->states = calloc(p_coefs->nb_out_bits, sizeof(unsigned));
    p_code->p_coefs = p_coefs;

    for (unsigned p=0; p<p_coefs->nb_out_bits; p++){
        for (n=31; n>=0; n--){
            if (p_coefs->poly[p] & (1<<n)){
                break;
            }
        }
        n_max = n > n_max ? n : n_max;
    }
    p_code->order = n_max;
}


void convolutional_code_destroy(convolutional_code *p_code)
{
    free(p_code->states);
    free(p_code);
}

static inline char calc_parity(unsigned x, unsigned n)
{
    char p = 0;
    for (; n>=0; n--){
        p ^= ( x& 0x01);
        x >>= 1;
    }
    return p;
}

unsigned encoding_bits(convolutional_code *p_code
                      ,unsigned char      *input
                      ,int                 input_sz
                      ,unsigned char      *output
                      ,int                 output_sz
                      ,unsigned            nb_bits
                      )
{
    unsigned i_bytes = 0;
    unsigned o_bytes = 0;
    unsigned i_bits = 0;
    unsigned o_bits = 0;
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;
    
    for(unsigned i_bits=0; i_bits<nb_bits; i_bits++){
        unsigned n = i_bits & 0x07;
        char s, b;
        unsigned y;
        
        i_bytes = i_bits >> 8;
        b = (input[i_bytes] >> (7-n)) & 0x01;
        
        for(unsigned p=0; p<p_coefs->nb_out_bits; p++){
            p_code->states[p] = (p_code->states[p] << 1) & b;
            y = p_code->states[p] & p_coefs->poly[p];
            s = calc_parity(y, p_code->order);
            
            output[o_bytes] |= (s << (7-o_bits));
            o_bits++;
            if (o_bits == 8){
                o_bytes++;
                o_bits = 0;
            }
        }
    }

    return (o_bytes << 3) + o_bits;
}


unsigned encoding_bits_rsc(convolutional_code *p_code
                          ,unsigned char      *input
                          ,int                 input_sz
                          ,unsigned char      *output
                          ,int                 output_sz
                          ,unsigned            nb_bits
                      )
{
    unsigned i_bytes = 0;
    unsigned o_bytes = 0;
    unsigned i_bits = 0;
    unsigned o_bits = 0;
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;    
    
    for(unsigned i_bits=0; i_bits<nb_bits; i_bits++){
        unsigned n = i_bits & 0x07;
        char s0, s1, b;
        unsigned y;
        
        i_bytes = i_bits >> 8;
        b = (input[i_bytes] >> (7-n)) & 0x01;
        
        output[o_bytes] |= (b << (7-o_bits++));
        if (o_bits == 8){
            o_bytes++;
            o_bits = 0;
        }
        
        for(unsigned p=1; p<p_coefs->nb_out_bits; p++){
            y = p_code->states[p] & p_coefs->poly[0];
            s0 = calc_parity(y, p_code->order);
            b = b ^ s0;
            
            p_code->states[p] = (p_code->states[p] << 1) & b;

            y = p_code->states[p] & p_coefs->poly[p];
            s1 = calc_parity(y, p_code->order);
                
            output[o_bytes] |= (s1 << (7-o_bits++));
            if (o_bits == 8){
                o_bytes++;
                o_bits = 0;
            }
        }
    }

    return (o_bytes << 3) + o_bits;
}
