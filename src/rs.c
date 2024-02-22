/* Reed-Soloman codes
 * here we implement punctualed RS code, remove 0 from the locator.
 */


#include "gf2m.h"
#include "gf2.h"
#include "rs.h"

struct rs_code_s{
    unsigned k;
    unsigned n;
    unsigned d;
    
    uint8_t **generator;
    uint8_t *genpoly;
    uint8_t *dec_mat_0;
    uint8_t *dec_mat_1;
    const uint8_t *alphas;

    uint8_t *scratch;
    
};

/* compute the generator polynomials */
static void rs_code_genpoly(rs_code *enc)
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


rs_code* rs_code_construct(int n, int k)
{
    rs_code *p_state = calloc(sizeof(rs_code), 1);
    int m = deg_plus_one(n);

    
    p_state->n = n;
    p_state->k = k;
    p_state->d = n - k + 1;
    
    assert((n-k) % 2 == 0);
    
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
    rs_code_genpoly(p_state);


    /* the decoding matrix, the first k+t columns are fixed like Vandermont
     * matrix, the following t+1 columns are fixed power of the primitive
     * elements times the recved words. here t=(n-k)/2 and shall be even
     */
    {
	int t = (n-k)/2;
	p_state->dec_mat_0 = (uint8_t*)calloc(sizeof(uint8_t), n*n);
	p_state->dec_mat_1 = (uint8_t*)calloc(sizeof(uint8_t), n*n);
	
	for (unsigned i=0; i<n; i++){
	    unsigned j=0; 
	    for(; j<k+t; j++){
		p_state->dec_mat_0[i*n + j] = p_state->dec_mat_1[i*n + j]  = gf2m_power_nz(p_state->alphas[i], j);
	    }
	    for(unsigned p=0; p<t; p++, j++){
		p_state->dec_mat_0[i*n + j] = gf2m_power_nz(p_state->alphas[i], p+1);
	    }
	}
    }

    p_state->scratch = malloc(4*n);
    
    return p_state;
}

int rs_get_genpoly(rs_code *enc, uint8_t* out, int out_sz)
{
    assert(out_sz >= enc->d);
    memcpy(out, enc->genpoly, enc->d);
    return enc->d;
}

void rs_code_destroy(rs_code *p)
{
    if (p){
        for (unsigned i=0; i<p->k; i++){
            free(p->generator[i]);
        }
        free(p->generator);
	free(p->genpoly);
	free(p->dec_mat_0);
	free(p->dec_mat_1);
	free(p->scratch);
    }
    free(p);
}


unsigned rs_encode(rs_code *p_enc
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

unsigned rs_encode_sys(rs_code *p_enc
		      ,unsigned char      *in
		      ,int              in_sz
		      ,unsigned char      *out
		      ,int              out_sz
		      )
{
    const unsigned n  = p_enc->n;
    const unsigned k  = p_enc->k;
    int dq, dr;
    uint8_t *in1 = p_enc->scratch;
    uint8_t *q = p_enc->scratch + n;
    memcpy(&in1[n-k], in, k);

    gf2m_poly_div(in1, n, p_enc->genpoly, n-k+1, q, n, &dq, &dr);
    assert(dr == n-k);

    memcpy(out, in1, n-k);
    memcpy(&out[n-k], in ,k);

    return n;
}



unsigned rs_decode_ge(rs_code         *p_rs
		     ,unsigned char   *in
		     ,int              in_sz
		     ,unsigned char   *out
		     ,int              out_sz
		     )
{

    const unsigned t = (p_rs->n - p_rs->k) / 2;
    const unsigned n  = p_rs->n;
    const unsigned k  = p_rs->k;
    int dq = 1, dr = 1;
    uint8_t *b = p_rs->scratch;
    uint8_t *lambda = b + n;
    
    for (unsigned i=0; i<n; i++){
	unsigned j =0;
	for (j=0; j<k+t; j++){
	    p_rs->dec_mat_1[i*n + j] = p_rs->dec_mat_0[i*n + j];
	}
	for (; j<n; j++){
	    p_rs->dec_mat_1[i*n + j] = gf2m_mult(in[i], p_rs->dec_mat_0[i*n + j]);
	}
	b[i] = in[i];
    }

    gf2m_gaussian_elemination(p_rs->dec_mat_1, n, n, b, n);

    /* now b shall contain the solution */
    lambda[0] = 1;
    memcpy(&lambda[1], &b[k+t], t);

    gf2m_poly_div(b, k+t, lambda, t+1, out, out_sz, &dq, &dr);

    if (dr == 0){
	/* lambda(x) is a divisor, */
	assert(dq == k);
	return k;
    }

    /* decoding failed */
    return 0;
}
