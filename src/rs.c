/* Reed-Soloman codes
 * here we implement punctualed RS code, remove 0 from the locator.
 */


#include "gf2_8.h"
#include "rs.h"

struct rs255_encoder_s{
    unsigned k;
    uint8_t **generator;
};

const unsigned rs255_n = 255;

rs255_encoder* rs255_encoder_construct(int k)
{
    rs255_encoder *p_state = calloc(sizeof(rs255_encoder), 1);
    p_state->k = k;
    p_state->generator = (uint8_t**)calloc(sizeof(uint8_t*), k);

    for (unsigned i=0; i<k; i++){
        uint8_t beta = alphas[i];
        p_state->generator[i] = (uint8_t*)calloc(sizeof(uint8_t), rs255_n);

        for (unsigned j=0; j<rs255_n; j++){
            p_state->generator[i][j] = gf2m_power_nz(beta, j);
        }
    }
    return p_state;
}

void rs255_encoder_destroy(rs255_encoder *p)
{
    if (p){
        for (unsigned i=0; i<p->k; i++){
            free(p->generator[i]);
        }
        free(p->generator);
    }
    free(p);
}

unsigned rs255_encode(rs255_encoder *p_enc
                     ,unsigned char      *in
                     ,int              in_sz
                     ,unsigned char      *out
                     ,int              out_sz
                     )
{
    assert(in_sz == p_enc->k);
    assert(out_sz >= rs255_n);

    (void)in_sz;
    (void)out_sz;

    for (unsigned i = 0; i<rs255_n; i++){
        uint8_t fi = 0;
        for (unsigned j=0; j<p_enc->k; j++){
            fi = gf2m_add(fi, gf2m_mult(in[j], p_enc->generator[j][i]));
        }
        out[i] = fi;
    }
    return rs255_n;
}

