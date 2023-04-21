#include "conv.h"
#include <string.h>
#include <stdint.h>
#include <assert.h>

typedef  struct{
    uint16_t curr_s;
    uint16_t prev_s[2]; /* previous state if input is 0 or 1 */
    char     obits[2]; /* maximum 16 bits output  */
}trellis_node;


struct convolutional_code_s
{
    convolutional_code_coefs *p_coefs;
    unsigned state;
    unsigned order; //K -1
    trellis_node *trellis;

    /* viterbi decoding */
    float *p_metrics; //current running metrics, each state has the best to date metric
    unsigned delay_len;
    unsigned delay_idx;
    uint16_t **pp_path; // history of survivors.
    char **pp_b_active; // whether the state is active
    char *p_b_merged;   //whehter the path has merged
    uint16_t *p_active_state;
};

#define OUTPUT_ONE_BIT(out, obyte, obit, bit) do{     \
        out[obyte] |= ((bit) << (7-(obit++)));    \
        if (obit == 8){                              \
            obit = 0;                                 \
            out[++obyte] = 0;                         \
        }                                             \
    }while(0)

static inline char calc_parity(unsigned x, unsigned n)
{
    char p = 0;
    for (unsigned i=0; i<=n; i++){
        p ^= ( x& 0x01);
        x >>= 1;
    }
    return p;
}
    
convolutional_code* convolutional_code_construct( convolutional_code_coefs *p_coefs)
{
    convolutional_code *p_code = calloc(1, sizeof(convolutional_code));
    unsigned n_max = 0, msb, max_s, n;

    /* rate 1/7 is hte lowest*/
    assert(p_coefs->nb_out_bits < 8);
    
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
    assert(n_max < 32); /* maixmum contraints lenght 32 */

    if (n_max < 12) {
        /* when constraint length is less or equal than 12, K<=12, 
         * we can enable viterbi decoding */
        
        /* build the trellis, the output bits is fully determined by (s[k], s[k+1]) 
         * e.g. s[k+1] = 10010, then previous state could be 11001 or 01001, but the 
         * input can only be zero.  
         */
        msb = 1 << (p_code->order-1);
        max_s = 1<<(p_code->order);
        p_code->trellis = calloc(max_s, sizeof(trellis_node));
        for (unsigned i=0; i<max_s; i++){
            const char b = i & 0x01;
            unsigned s = i >> 1;
            p_code->trellis[i].curr_s = i;
            p_code->trellis[i].prev_s[0] = s;
            p_code->trellis[i].prev_s[1] = s | msb;
            if (p_coefs->b_rsc){
                const unsigned fb_tap = (p_coefs->poly[0] & (~1));
                for (unsigned j = 0; j<2; j++){
                    s = (p_code->trellis[i].prev_s[j] << 1) | b;
                    p_code->trellis[i].obits[j] = b ^ calc_parity(s & fb_tap, p_code->order);
                    for(unsigned p=1; p<p_coefs->nb_out_bits; p++){
                        unsigned y = s & p_coefs->poly[p];
                        p_code->trellis[i].obits[j] = (p_code->trellis[i].obits[j] << 1) | calc_parity(y, p_code->order);
                    }
                }
            }
            else{
                for (unsigned j = 0; j<2; j++){
                    /* this is input + state */
                    s = (p_code->trellis[i].prev_s[j] << 1) | b;
                    for(unsigned p=0; p<p_coefs->nb_out_bits; p++){
                        unsigned y = s & p_coefs->poly[p];
                        p_code->trellis[i].obits[j] = (p_code->trellis[i].obits[j] << 1) | calc_parity(y, p_code->order);
                    }
                }
            }
        }

        p_code->p_metrics = calloc(sizeof(float), max_s);
        p_code->delay_len = p_coefs->dec_delay_factor * (p_code->order+1);
        p_code->pp_path = calloc(sizeof(void*), p_code->delay_len);
        p_code->pp_active = calloc(sizeof(void*), p_code->delay_len);
        for (unsigned j=0; j<p_code->delay_len; j++){
            p_code->pp_path[j] = calloc(sizeof(uint16_t), max_s);
            p_code->pp_active[j] = calloc(sizeof(char), max_s);
        }
        p_code->p_b_merged = calloc(sizeof(char), p_code->delay_len);
        p_code->p_active_state = calloc(sizeof(uint16_t), p_code->delay_len);
    }
    return p_code;
}


void convolutional_code_destroy(convolutional_code *p_code)
{
    free(p_code->trellis);
    if (p_code->pp_path){
        for (unsigned j=0; j<p_code->delay_len; j++){
            free(p_code->pp_path[j]);
        }
        free(p_code->pp_path);        
    }

    if (p_code->pp_active){
        for (unsigned j=0; j<p_code->delay_len; j++){
            free(p_code->pp_active[j]);
        }
        free(p_code->pp_active);        
    }

    free(p_code->p_b_merged);
    free(p_code->p_active_state);
    free(p_code->p_metrics);
    free(p_code);
}


void get_trellis_node(convolutional_code *p_code
                     ,unsigned node
                     ,unsigned *p0
                     ,unsigned *p1
                     ,unsigned      *p_in_bit
                     ,unsigned      *p_out0_bits
                     ,unsigned      *p_out1_bits
                     )

{
    if(p_code->trellis && p0 && p1 && p_in_bit && p_out0_bits && p_out1_bits)
    {
        *p0 = p_code->trellis[node].prev_s[0];
        *p1 = p_code->trellis[node].prev_s[1];
        *p_in_bit = node & 0x01;
        *p_out0_bits = p_code->trellis[node].obits[0];
        *p_out1_bits = p_code->trellis[node].obits[1];
    }
}


static unsigned encode_bits(convolutional_code *p_code
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

    output[0] = 0;
    for(unsigned i_bits=0; i_bits<nb_bits; i_bits++){
        unsigned n = i_bits & 0x07;
        char s, b;
        unsigned y;
        
        i_bytes = i_bits >> 3;
        b = (input[i_bytes] >> (7-n)) & 0x01;
        p_code->state = (p_code->state << 1) | b;        
        
        for(unsigned p=0; p<p_coefs->nb_out_bits; p++){
            y = p_code->state & p_coefs->poly[p];
            s = calc_parity(y, p_code->order);
            OUTPUT_ONE_BIT(output, o_bytes, o_bits, s);
        }
    }
    /* terminal the code */
    for(unsigned i_bits=0; i_bits<=p_code->order; i_bits++){
        p_code->state = p_code->state << 1;
        char s;
        unsigned y;
        for(unsigned p=0; p<p_coefs->nb_out_bits; p++){
            y = p_code->state & p_coefs->poly[p];
            s = calc_parity(y, p_code->order);
            OUTPUT_ONE_BIT(output, o_bytes, o_bits, s);            
        }
    }

    return (o_bytes << 3) + o_bits;
}


static unsigned encode_bits_rsc(convolutional_code *p_code
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
    const unsigned fb_tap = (p_coefs->poly[0] & (~1));

    output[0] = 0;    
    for(unsigned i_bits=0; i_bits<nb_bits; i_bits++){
        unsigned n = i_bits & 0x07;
        char s0, s1, b;
        unsigned y;
        
        i_bytes = i_bits >> 3;
        b = (input[i_bytes] >> (7-n)) & 0x01;
        
        OUTPUT_ONE_BIT(output, o_bytes, o_bits, b);
        
        y = p_code->state & fb_tap;
        s0 = calc_parity(y, p_code->order);
        b = b ^ s0;
        p_code->state = (p_code->state<< 1) | b;
        
        for(unsigned p=1; p<p_coefs->nb_out_bits; p++){
            y = p_code->state & p_coefs->poly[p];
            s1 = calc_parity(y, p_code->order);
            OUTPUT_ONE_BIT(output, o_bytes, o_bits, s1);
        }
    }

    /* terminating the code */
    for(unsigned i_bits=0; i_bits<=p_code->order; i_bits++){
        char s0,s1;
        unsigned y;
        y = p_code->state & fb_tap;
        s0 = calc_parity(y, p_code->order);
        OUTPUT_ONE_BIT(output, o_bytes, o_bits, s0);
        
        p_code->state = (p_code->state<< 1);
        for(unsigned p=1; p<p_coefs->nb_out_bits; p++){
            y = p_code->state & p_coefs->poly[p];
            s1 = calc_parity(y, p_code->order);
            OUTPUT_ONE_BIT(output, o_bytes, o_bits, s1);            

        }
    }
    return (o_bytes << 3) + o_bits;
}


unsigned encode(convolutional_code *p_code
               ,unsigned char      *in
               ,int              in_sz
               ,unsigned char      *out
               ,int              out_sz
               ,unsigned            nb_bits
               )
{
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;
    if (p_coefs->b_rsc){
        return encode_bits_rsc(p_code, in, in_sz, out, out_sz, nb_bits);
    }
    return encode_bits(p_code, in, in_sz, out, out_sz, nb_bits);
}


/* r[0] r[1] r[2] r[3]  corresponds to b3b2b1b0 */
static float inner_prod(float *in, char bits, unsigned nb)
{
    float sum = 0.0f;
    for (unsigned i=0; i<nb; i++){
        float a = ((bits >> (nb-i-1)) & 0x1) - 1.0f;
        sum += in[i] * a;
    }
}

unsigned viterbi_decode(convolutional_code *p_code
                       ,float              *input
                       ,int                 input_sz
                       ,unsigned char      *output
                       ,int                 output_sz
                       ,unsigned            nb_bits
                      )
{
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;
    unsigned max_s = 1<<p_code->order;
    const unsigned idx0 = p_code->delay_idx;
    unsigned obytes = 0, obits = 0;
    if (!p_code->trellis){
        /* K is too large to use viterbi decoding */
        return 0 ;
    }
    output[0] = 0;

    for (unsigned i = 0;  i<nb_bits; i+= p_coefs->nb_out_bits){

        //go through all the states
        for (unsigned s = 0; s<max_s; s++){

            float m0 = inner_prod(&input[i], p_code->trellis[s].obits[0], p_coefs->nb_out_bits);
            float m1 = inner_prod(&input[i], p_code->trellis[s].obits[1], p_coefs->nb_out_bits);

            p_code->pp_path[idx0][s] = m0 > m1 ? p_code->trellis[s].prev_s[0] : p_code->trellis[s].prev_s[1];
            p_code->pp_b_active[idx0][s] = 1;
            p_code->p_metrics[s] += m0 > m1 ? m0 : m1;
        }
        p_code->p_b_merged[idx0] = 1;
            
        //start to output decoded bits 
        if (i >= (p_code->delay_len - 1 )*p_coefs->nb_out_bits){
            unsigned oldest = (idx0 + 1) % p_code->delay_len;
            
            if (!p_code->p_b_merged[oldest]){
                for (unsigned k =0; k<p_code->delay_len; k++){
                    unsigned idx  = (idx0 - k + p_code->delay_len) % p_code->delay_len;
                    unsigned idx1 = (idx- 1 + p_code->delay_len) % p_code->delay_len;
                    char b_merged = 1;
                    char b_init = 0;
                
                    for (unsigned s = 0; s<max_s; s++){
                        if (p_code->pp_b_active[idx][s]){
                            uint16_t prev_s = p_code->pp_path[idx][s];
                            p_code->pp_b_active[idx1][prev_s] = 1;
                            if (!b_init){
                                b_init = 1;
                            }
                            else if (p_code->p_active_state[idx1] != prev_s){
                                b_merged = 0;
                            }
                            p_code->p_active_state[idx1] = prev_s;                            
                        }
                    }
                    p_code->p_b_merged[idx1] = b_merged;
                }
            }

            if (!p_code->p_b_merged[oldest]){
                //still not converged, just randomly pick
                printf("%s \n", "path nor merged");
            }
            else{
                uint16_t curr_s = p_code->p_active_state[oldest];
                char b;
                if (p_coefs->b_rsc){
                    uint16_t prev_s = p_code->pp_path[oldest][curr_s];
                    char msb = (prev_s >> (p_code->order -1 )) & 0x01 ;
                    b = (curr_s & 0x01) ^ msb;
                }
                else{
                    b = curr_s & 0x01;
                }
                OUTPUT_ONE_BIT(output, obytes, obits, b);
            }
        }
    }

    /*epilog, as we know the last state is 0 state, we trace back
      from this state. */
    for (unsigned i=0; i<p_code->delay_len; i++){
        


    }
        
        
    
    



}
