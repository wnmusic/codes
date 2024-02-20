#include "conv.h"
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "debug_trace.h"
#include <math.h>

typedef  struct{
    int16_t curr_s;
    int16_t prev_s[2]; /* previous state if input is 0 or 1 */
    char    obits[2]; /* maximum 16 bits output  */
}trellis_node;

struct convolutional_code_s
{
    convolutional_code_coefs *p_coefs;
    unsigned state;
    unsigned order; //K -1
    trellis_node *trellis;

    /* viterbi decoding */
    float *p_metrics; //current running metrics, each state has the best to date metric
    float *p_prev_metrics; //previous metrics
    
    unsigned delay_len;
    int16_t **pp_path; // history of survivors.
    char **pp_b_active; // whether the state is active
    char *p_b_merged;   //whehter the path has merged
    int16_t *p_active_state;


    /* fano deocding */
    float curr_snr; /* currnet estimated signal to noise ratio, might change between blocks  */
    float *p_metrics_table;
    float th_delta;    
    float running_metric;
    float running_th;
    int   new_flag;
    int   bit_index;
    
};

typedef struct
{
    unsigned char *buf;
    unsigned obyte;
    unsigned obit;
}bitstream;
    
#define OUTPUT_ONE_BIT(out, bit) do{                            \
        out->buf[out->obyte] |= ((bit) << (out->obit++));       \
        if (out->obit == 8){                                    \
            out->obit = 0;                                      \
            out->buf[++out->obyte] = 0;                         \
        }                                                       \
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
                        p_code->trellis[i].obits[j] = p_code->trellis[i].obits[j] | (calc_parity(y, p_code->order) << p);
                    }
                }
            }
            else{
                for (unsigned j = 0; j<2; j++){
                    /* this is input + state */
                    s = (p_code->trellis[i].prev_s[j] << 1) | b;
                    for(unsigned p=0; p<p_coefs->nb_out_bits; p++){
                        unsigned y = s & p_coefs->poly[p];
                        p_code->trellis[i].obits[j] = p_code->trellis[i].obits[j] | (calc_parity(y, p_code->order) << p);
                    }
                }
            }
        }

        p_code->p_metrics = calloc(sizeof(float), max_s);
        p_code->p_prev_metrics = calloc(sizeof(float), max_s);
        p_code->delay_len = p_coefs->dec_delay_factor * (p_code->order+1);
        p_code->pp_path = calloc(sizeof(void*), p_code->delay_len);
        p_code->pp_b_active = calloc(sizeof(void*), p_code->delay_len);
        for (unsigned j=0; j<p_code->delay_len; j++){
            p_code->pp_path[j] = calloc(sizeof(int16_t), max_s);
            p_code->pp_b_active[j] = calloc(sizeof(char), max_s);
        }
        p_code->p_b_merged = calloc(sizeof(char), p_code->delay_len);
        p_code->p_active_state = calloc(sizeof(int16_t), p_code->delay_len);
        for (unsigned i=0; i<p_code->delay_len; i++){
            p_code->p_active_state[i] = -1;
        }
    }

    DBG_TRACE_OPEN();
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

    if (p_code->pp_b_active){
        for (unsigned j=0; j<p_code->delay_len; j++){
            free(p_code->pp_b_active[j]);
        }
        free(p_code->pp_b_active);        
    }

    free(p_code->p_b_merged);
    free(p_code->p_active_state);
    free(p_code->p_metrics);
    free(p_code->p_prev_metrics);
    free(p_code);

    DBG_TRACE_CLOSE();
}


void get_trellis_node(convolutional_code *p_code
                     ,unsigned            node
                     ,unsigned            *p0
                     ,unsigned            *p1
                     ,unsigned            *p_in_bit
                     ,unsigned            *p_out0_bits
                     ,unsigned            *p_out1_bits
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


static void encode_bits(convolutional_code *p_code
                       ,unsigned char      *input
                       ,bitstream          *output
                       ,unsigned            nb_bits
                       )
{
    unsigned i_bytes = 0;
    unsigned i_bits = 0;
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;

    for(unsigned i_bits=0; i_bits<nb_bits; i_bits++){
        unsigned n = i_bits & 0x07;
        char s, b;
        unsigned y;
        
        i_bytes = i_bits >> 3;
        b = (input[i_bytes] >> n) & 0x01;
        p_code->state = (p_code->state << 1) | b;        
        
        for(unsigned p=0; p<p_coefs->nb_out_bits; p++){
            y = p_code->state & p_coefs->poly[p];
            s = calc_parity(y, p_code->order);
            OUTPUT_ONE_BIT(output, s);
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
            OUTPUT_ONE_BIT(output, s);            
        }
    }

}


static void encode_bits_rsc(convolutional_code *p_code
                               ,unsigned char      *input
                               ,bitstream          *output
                               ,unsigned            nb_bits                           
                               )
{
    unsigned i_bytes = 0;
    unsigned i_bits = 0;
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;
    const unsigned fb_tap = (p_coefs->poly[0] & (~1));

    for(unsigned i_bits=0; i_bits<nb_bits; i_bits++){
        unsigned n = i_bits & 0x07;
        char s0, s1, b;
        unsigned y;
        
        i_bytes = i_bits >> 3;
        b = (input[i_bytes] >> n) & 0x01;
        
        OUTPUT_ONE_BIT(output,  b);
        
        y = (p_code->state << 1) & fb_tap;
        s0 = calc_parity(y, p_code->order);
        b = b ^ s0;
        p_code->state = (p_code->state<< 1) | b;
        
        for(unsigned p=1; p<p_coefs->nb_out_bits; p++){
            y = p_code->state & p_coefs->poly[p];
            s1 = calc_parity(y, p_code->order);
            OUTPUT_ONE_BIT(output, s1);
        }
    }

    /* terminating the code */
    for(unsigned i_bits=0; i_bits<=p_code->order; i_bits++){
        char s0,s1;
        unsigned y;
        y = (p_code->state<<1) & fb_tap;
        s0 = calc_parity(y, p_code->order);
        OUTPUT_ONE_BIT(output, s0);
        
        p_code->state = (p_code->state<< 1);
        for(unsigned p=1; p<p_coefs->nb_out_bits; p++){
            y = p_code->state & p_coefs->poly[p];
            s1 = calc_parity(y, p_code->order);
            OUTPUT_ONE_BIT(output, s1);            
        }
    }
}


unsigned encode(convolutional_code *p_code
               ,unsigned char      *in
               ,int                 in_sz
               ,unsigned char      *out
               ,int                 out_sz
               ,unsigned            nb_bits
               )
{
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;
    bitstream stream;
    stream.buf = out;
    stream.obyte = 0;
    stream.obit = 0;
    stream.buf[0] = 0;
    
    (void)in_sz;
    (void)out_sz;
    
    if (p_coefs->b_rsc){
        encode_bits_rsc(p_code, in, &stream, nb_bits);
    }else{
        encode_bits(p_code, in, &stream, nb_bits);
    }

    return stream.obyte*8 + stream.obit;
}


/* for bits inside a byte, the mapping order will be 
 * r[0] r[1] ... r[7] maps to bits b0, b1, ..., b7 of a byte 
 */
static float bpsk_inner_prod(float *in, char bits, unsigned nb)
{
    float sum = 0.0f;
    for (unsigned i=0; i<nb; i++){
        float a = 1.0f - 2.0f *((bits >> i) & 0x1);
        sum += in[i] * a;
    }
    return sum;
}

static inline
void 
decode_one_bit(convolutional_code *p_code
              ,unsigned output_idx
              ,bitstream *output
              )
{
    int16_t curr_s = p_code->p_active_state[output_idx];
    char b;
    if (p_code->p_coefs->b_rsc){
        int16_t prev_s = p_code->pp_path[output_idx][curr_s];
        const unsigned fb_tap = (p_code->p_coefs->poly[0] & (~1));        
        char msb = calc_parity( (prev_s << 1 ) & fb_tap, p_code->order);
        b = (curr_s & 0x01) ^ msb;
    }
    else{
        b = curr_s & 0x01;
    }
    OUTPUT_ONE_BIT(output, b);
}

unsigned viterbi_decode(convolutional_code *p_code
                       ,float              *input
                       ,int                 input_sz
                       ,unsigned char      *out
                       ,int                 out_sz
                       ,unsigned            nb_bits    
                      )
{
    const convolutional_code_coefs *p_coefs = p_code->p_coefs;
    unsigned max_s = 1<<p_code->order;
    unsigned idx0 = 0;
    unsigned obytes = 0, obits = 0;
    unsigned output_idx = 0;
    bitstream stream;
    int16_t curr_s = 0;
    
    (void)input_sz;
    (void)out_sz;

    if (!p_code->trellis){
        /* K is too large to use viterbi decoding */
        return 0 ;
    }
    out[0]       = 0;
    stream.buf   = out;
    stream.obyte = 0;
    stream.obit  = 0;

    for (unsigned i = 0;  i<nb_bits; i+= p_coefs->nb_out_bits){

        //go through all the states
        float min_m = 1e12f;
        for (unsigned s = 0; s<max_s; s++){

            float m0 = bpsk_inner_prod(&input[i], p_code->trellis[s].obits[0], p_coefs->nb_out_bits);
            float m1 = bpsk_inner_prod(&input[i], p_code->trellis[s].obits[1], p_coefs->nb_out_bits);
            int16_t prev_s0 = p_code->trellis[s].prev_s[0];
            int16_t prev_s1 = p_code->trellis[s].prev_s[1];

            m0 += p_code->p_prev_metrics[prev_s0];
            m1 += p_code->p_prev_metrics[prev_s1];

            p_code->pp_path[idx0][s] = m0 > m1 ?  prev_s0 :  prev_s1;
            p_code->pp_b_active[idx0][s] = 1;
            p_code->p_metrics[s] = m0 > m1 ? m0 : m1;
            
            min_m = p_code->p_metrics[s] < min_m ? p_code->p_metrics[s] : min_m;
        }

        //normalize
        if (min_m  > 200.0f){
            for (unsigned s = 0; s<max_s; s++){
                p_code->p_prev_metrics[s] = p_code->p_metrics[s] - min_m;
            }
        }
        else{
            for (unsigned s = 0; s<max_s; s++){
                p_code->p_prev_metrics[s] = p_code->p_metrics[s];
            }
        }            

        p_code->p_b_merged[idx0] = 0;

        DBG_TRACE_ARRAY_1D("metric", float, p_code->p_metrics, max_s);
        DBG_TRACE_ARRAY_1D("path", int16_t, p_code->pp_path[idx0], max_s);
            
        //start to output decoded bits 
        if (i >= (p_code->delay_len - 1 )*p_coefs->nb_out_bits){
            unsigned idx, idx1;
            if (!p_code->p_b_merged[output_idx]){
                /* as we are determine if previous state has merged, the loop
                   should end when previous state is the oldest */
                for (unsigned k =0; k<p_code->delay_len - 1; k++){
                    idx  = (idx0 - k + p_code->delay_len) % p_code->delay_len;
                    idx1 = (idx- 1 + p_code->delay_len) % p_code->delay_len;
                    
                    if (p_code->p_b_merged[idx]){
                        curr_s = p_code->p_active_state[idx];
                        p_code->p_b_merged[idx1] = 1;
                        p_code->p_active_state[idx1] = p_code->pp_path[idx][curr_s];   
                    }
                    else{
                        char b_merged = 1;
                        memset(p_code->pp_b_active[idx1], 0, sizeof(char)*max_s);
                        for (unsigned s = 0; s<max_s; s++){
                            if (p_code->pp_b_active[idx][s]){
                                int16_t prev_s = p_code->pp_path[idx][s];
                                p_code->pp_b_active[idx1][prev_s] = 1;
                                if (p_code->p_active_state[idx1] != -1 && p_code->p_active_state[idx1] != prev_s){
                                    b_merged = 0;
                                }
                                p_code->p_active_state[idx1] = prev_s;                            
                            }
                        }
                        p_code->p_b_merged[idx1] = b_merged;
                    }
                }
            }

            if (!p_code->p_b_merged[output_idx]){
                //still not converged, just randomly pick, this will be a decoding error 
                //printf("output_idx: %d, idx %d, idx1: %d, i: %d, not merge \n", output_idx, idx, idx1, i);
            }

            decode_one_bit(p_code, output_idx, &stream);
            
            p_code->p_b_merged[output_idx] = 0;
            p_code->p_active_state[output_idx] = -1;

            output_idx = (output_idx + 1) % p_code->delay_len;
        }
        idx0 = (idx0 + 1) % p_code->delay_len;
    }

    /* we know that terminated state is 0 */
    for (unsigned i=0, curr_s = 0; i<p_code->delay_len; i++)
    {
        unsigned idx  = (idx0 - i - 1 + p_code->delay_len) % p_code->delay_len;
        if (p_code->p_b_merged[idx]){
            /* if the merged state is not what we traced back from ending state
             * a decoding error occurs. 

             */
            break;
        }
        p_code->p_b_merged[idx] = 1;
        p_code->p_active_state[idx] = curr_s;
        curr_s = p_code->pp_path[idx][curr_s];

        if (idx == output_idx){
            break;
        }
    }
    
    while(output_idx != idx0)
    {
        decode_one_bit(p_code, output_idx, &stream);
        output_idx = (output_idx + 1) % p_code->delay_len;
    }
        
    return stream.obyte*8 + stream.obit;
}        


unsigned fano_decode(convolutional_code *p_code
                    ,float              *input
                    ,int                 input_sz
                    ,unsigned char      *out
                    ,int                 out_sz
                    ,unsigned            nb_bits
                    ,float               snr_dB
                    )
{

    if (fabs(p_code->curr_snr - snr_dB) > 1.5f){
        p_code->curr_snr = snr_dB;
    }
    return 0;
}
