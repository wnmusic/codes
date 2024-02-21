/* Reed-Soloman codes
 * here we implement punctualed RS code, remove 0 from the locator.
 */


#include "gf2m.h"
#include "gf2.h"
#include "rs.h"

struct rs_encoder_s{
    unsigned k;
    unsigned n;
    unsigned d;
    
    uint8_t **generator;
    uint8_t *genpoly;
    const uint8_t *alphas;
    
};

/* compute the generator polynomials */
static void rs_encoder_genpoly(rs_encoder *enc)
{
    int d = enc->n - enc->k;
    uint8_t poly[2];
    uint8_t tmp[256];
    int curr_d; 
    memset(enc->genpoly, 0, d+1);
    enc->genpoly[0] = enc->alphas[1];
    enc->genpoly[1] = 1;
    curr_d = 2;
    for (int i=2; i<=d; i++){
        poly[0] = enc->alphas[i];
        poly[1] = 1;
        curr_d = gf2m_poly_mult(enc->genpoly, curr_d, poly, 2, tmp, curr_d+1);

        memcpy(enc->genpoly, tmp, curr_d);
    }
    return;
}


rs_encoder* rs_encoder_construct(int n, int k)
{
    rs_encoder *p_state = calloc(sizeof(rs_encoder), 1);
    int m = deg_plus_one(n);
    
    p_state->n = n;
    p_state->k = k;
    p_state->d = n - k + 1;
    
    p_state->alphas = setup_gf2m_ops(m);
    assert(p_state->alphas);
    
    p_state->generator = (uint8_t**)calloc(sizeof(uint8_t*), k);
    for (unsigned i=0; i<k; i++){
        uint8_t beta = p_state->alphas[i];
        p_state->generator[i] = (uint8_t*)calloc(sizeof(uint8_t), p_state->n);

        for (unsigned j=0; j<p_state->n; j++){
            p_state->generator[i][j] = gf2m_power_nz(beta, j);
        }
    }

    /* genpoly */
    p_state->genpoly = (uint8_t*)calloc(sizeof(uint8_t), p_state->d);
    rs_encoder_genpoly(p_state);
    return p_state;
}

int rs_get_genpoly(rs_encoder *enc, uint8_t* out, int out_sz)
{
    assert(out_sz >= enc->d);
    memcpy(out, enc->genpoly, enc->d);
    return enc->d;
}

void rs_encoder_destroy(rs_encoder *p)
{
    if (p){
        for (unsigned i=0; i<p->k; i++){
            free(p->generator[i]);
        }
        free(p->generator);
	free(p->genpoly);
    }
    free(p);
}


unsigned rs_encode(rs_encoder *p_enc
		  ,unsigned char      *in
		  ,int              in_sz
		  ,unsigned char      *out
		  ,int              out_sz
		  )
{
    assert(in_sz == p_enc->k);
    assert(out_sz >= p_enc->n);

    (void)in_sz;
    (void)out_sz;

    for (unsigned i = 0; i<p_enc->n; i++){
        uint8_t fi = 0;
        for (unsigned j=0; j<p_enc->k; j++){
            fi = gf2m_add(fi, gf2m_mult(in[j], p_enc->generator[j][i]));
        }
        out[i] = fi;
    }
    return p_enc->n;
}

