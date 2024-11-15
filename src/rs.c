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

    const uint8_t *alphas;
    const uint8_t *index_of_alphas;
    uint8_t  mul_order;

    int method;

    /* polynomrial evaluation construction */
    uint8_t **generator;
    uint8_t *dec_mat_0;
    uint8_t *dec_mat_1;
    
    /* BM method of decoding */
    uint8_t *genpoly;
    uint8_t *scratch;
    
};


/* compute the generator polynomials
 *  cyclic (n, k, d) RS code c(x) could be constructed by c(x) = g(x) m(x)
 *  where g(x) is a n-k degree polynomial m(x) is less than k. the degree of
 *  c(x) is less then n.
 *
 *  g(x) has zeors at a, a^2, a^(n-k) 
 */
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


rs_code* rs_code_construct(int n, int k, int method)
{
    rs_code *p_state = calloc(sizeof(rs_code), 1);
    int m = deg_plus_one(n);
    int z = 1;
    
    p_state->n = n;
    p_state->k = k;
    p_state->d = n - k + 1;
    
    assert((n-k) % 2 == 0);
    for (int i=0; i<m;i++){
	z *= 2;
    }
    assert( (n+1) == z && "make sure that n can be used as mask "); 
    
    setup_gf2m_ops(m, &p_state->alphas, &p_state->index_of_alphas, &p_state->mul_order);
    assert(p_state->alphas);

    if (method == RS_CODE_POLYVAL){
	
	p_state->generator = (uint8_t**)calloc(sizeof(uint8_t*), k);
	for (unsigned i=0; i<k; i++){
	    uint8_t beta = p_state->alphas[i];
	    p_state->generator[i] = (uint8_t*)calloc(sizeof(uint8_t), p_state->n);

	    for (unsigned j=0; j<p_state->n; j++){
		p_state->generator[i][j] = gf2m_power_nz(beta, j);
	    }
	}

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
    }
    else{
	
	/* genpoly */
	p_state->genpoly = (uint8_t*)calloc(sizeof(uint8_t), p_state->d);
	rs_code_genpoly(p_state);

	/* the decoding matrix, the first k+t columns are fixed like Vandermont
	 * matrix, the following t+1 columns are fixed power of the primitive
	 * elements times the recved words. here t=(n-k)/2 and shall be even
	 */
	p_state->scratch = malloc(8*n);
    }
    p_state->method = method;
    
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


static unsigned rs_encode_polyval(rs_code *p_enc
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

static unsigned rs_encode_sys(rs_code *p_enc
			     ,unsigned char      *in
			     ,int              in_sz
			     ,unsigned char      *out
			     ,int              out_sz
			     )
{
    const unsigned n  = p_enc->n;
    const unsigned k  = p_enc->k;
    uint8_t *g_idx = p_enc->scratch;

    /* the degree of generator polynomial is n-k, the number of parity check
       and it is monontic.
       here we use 255 to represent 0's power. 
     */
    memset(out, 0, n-k);
    memcpy(&out[n-k], in, k);
    
    for (int i = 0; i<n-k;i++){
	g_idx[i] = p_enc->index_of_alphas[ p_enc->genpoly[i] ];
    }
    
    for (int i=k-1; i>=0; i--){
	uint8_t m = in[i] ^ out[n-k-1];
	uint8_t m_idx = p_enc->index_of_alphas[m];

        for (int j=n-k-1; j>=1; j--){
	    if (m  && p_enc->genpoly[j]){
		out[j] = out[j-1] ^ p_enc->alphas[ (m_idx + g_idx[j] ) % p_enc->mul_order];
	    }else{
		out[j] = out[j-1];
	    }
        }
	if (m && p_enc->genpoly[0]){
	    out[0] = p_enc->alphas[ (m_idx + g_idx[0] ) % p_enc->mul_order];
	}else{
	    out[0] = 0;
	}
    }

    return n;
}

unsigned rs_encode(rs_code            *p_enc
		  ,unsigned char      *in
		  ,int                 in_sz
		  ,unsigned char      *out
		  ,int                 out_sz
		  )
{
    if (p_enc->method == RS_CODE_POLYVAL){
	return rs_encode_polyval(p_enc, in, in_sz, out, out_sz);
    }else if (p_enc->method == RS_CODE_CYCLIC){
	return rs_encode_sys(p_enc, in, in_sz, out, out_sz);
    }
    assert(0 && "not supported encoding method");
    return 0;
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

    if (p_rs->method != RS_CODE_POLYVAL){
	return 0;
    }
    
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
	return k;
    }

    /* decoding failed */
    return 0;
}

/* BM method to decode systematic RS codes */
    
/* calculate the synchrom
   in array stores the receiving symbols, the polynomnial will be
   r(x) = r0 + r_1X + r_2X^2 + + r_(n-1)
   the syndromes will just be  r(alpha), r(alpha^2) ... r(alpha^(n-k))
*/

#define EXPONET_ADD(m,n,order)  ((m) + (n)) >= (order) ? ((m)+(n)-(order)) : ((m)+(n))

static void rs_syndromes_calc(rs_code          *p_rs
			     ,unsigned char    *r
			     ,unsigned char    *s
			     ,uint8_t          *scratch
			     )
{
    const unsigned n  = p_rs->n;
    const unsigned k  = p_rs->k;
    const uint8_t order = p_rs->mul_order;

    uint8_t *r_idx = scratch;
    uint8_t *curr_idx = scratch + n;
    
    for (uint8_t i=0; i<n; i++){
	r_idx[i] = EXPONET_ADD(p_rs->index_of_alphas[ r[i] ], i, order);
    }

    for (uint8_t i=0; i<n-k; i++){
	uint8_t sum = 0;
	for (uint8_t j=0; j<n;j++){
	    if (r[j]){
		sum = sum ^ p_rs->alphas[ r_idx[j] ];
		r_idx[j] = EXPONET_ADD(r_idx[j], j, order);
	    }
	}
	s[i] = sum;
    }
}

/* this is Berlekamp-Massy interation to find the error locator poly
 * (1+b1X)*(1+b2X)*... *(1+btX)
 *  the input is syndrome sequence
 */
static void BM_iteration(rs_code       *p_rs
			,unsigned char *s        /* input syndrome */
			,int            nb_s			
			,unsigned char *lambda    /* output error locator polynomial*/
			,uint8_t       *scratch			
			)
{

    const unsigned n  = p_rs->n;
    const unsigned k  = p_rs->k;
    const unsigned t = (n - k) / 2;
    
    uint8_t *prev_lambda = scratch;
    uint8_t *tmp_lambda = scratch + t+1;    
    uint8_t L = 0;
    uint8_t x = 1;
    uint8_t b = 1;
    uint8_t tmp_L = 0;
    
    memset(lambda, 0, t+1);
    memset(prev_lambda, 0, t+1);
    lambda[0] = 1;
    prev_lambda[0] = 1;

    for (uint8_t i=0; i<nb_s; i++){
	uint8_t d = s[i];
	uint8_t scale;
	for (uint8_t j=1; j<=L && j<=i; j++){
	    d = d ^ gf2m_mult(lambda[j], s[i-j]);
	}
	if (d == 0){
	    x = x + 1;
	    continue;
	}

	scale = gf2m_div_nz(d,b);
	if (2*L > i){ /* we don't have to increase the degree of connection */

	    for (uint8_t j=x; j<=L;j++){
		lambda[j] = lambda[j] ^ gf2m_mult(scale, prev_lambda[j-x]);
	    }
	    x = x + 1;
	}else{
	    /* we now have to increase the degree, and shall store current
	       lambda to prev for later use
	    */
	    memcpy(tmp_lambda, lambda, L+1);
	    tmp_L = L;
	    L = i+1-L;
	    for (uint8_t j=x; j<=L;j++){
		lambda[j] = lambda[j] ^ gf2m_mult(scale, prev_lambda[j-x]);
	    }
	    memcpy(prev_lambda, tmp_lambda, tmp_L+1);
	    b = d;
	    x = 1;
	}
    }
}

unsigned rs_decode_BM(rs_code         *p_rs
		     ,unsigned char   *in
		     ,int              in_sz
		     ,unsigned char   *out
		     ,int              out_sz
		     )
{
    const unsigned t = (p_rs->n - p_rs->k) / 2;
    const unsigned n  = p_rs->n;
    const unsigned k  = p_rs->k;
    const uint8_t order = p_rs->mul_order;
    
    unsigned offset = 0;
    uint8_t	*syndroms      = p_rs->scratch + offset;  offset += n;   /* this will be reused, so n bytes */
    uint8_t	*err_poly      = p_rs->scratch + offset;  offset += t+1;
    uint8_t	*err_coefs_idx = p_rs->scratch + offset;  offset += t+1;
    uint8_t	*beta_idx      = p_rs->scratch + offset;  offset += t+1;
    uint8_t	*z_poly	       = p_rs->scratch + offset;  offset += t+1;
    uint8_t	*z_poly_idx    = p_rs->scratch + offset;  offset += t+1;
    uint8_t	*err_mag       = p_rs->scratch + offset;  offset += t+1;
    uint8_t p = 0;
    
    assert(offset < 8*n);
    
    rs_syndromes_calc(p_rs, in, syndroms, p_rs->scratch + offset);
    BM_iteration(p_rs, syndroms, n-k, err_poly, p_rs->scratch + offset);

    for (uint8_t i=0; i<=t; i++){
	err_coefs_idx[i] = p_rs->index_of_alphas[ err_poly[i] ];
    }

    for (int i=n-1; i>=0; i--){
	/* 1 test if i bit has errors, this is the same as testing
	 * if alpha^i is a root of err locator polynomial;
	 */
	uint8_t check = 1;
	for (uint8_t j=1; j<=t; j++){
	    if (err_poly[j]){
		err_coefs_idx[j] = EXPONET_ADD(err_coefs_idx[j], j, order);
		check = check ^ p_rs->alphas[ err_coefs_idx[j] ];
	    }
	}

	if (!check){
	    beta_idx[p]  = i;
	    p++;
	}
    }

    if (p > t){
	/* decoding failure */
	return 0;
    }
    else if (p==0){
	memcpy(out, &in[n-k], k);
	return k;
    }

    /* compute Z(x)  */
    z_poly[0] = 1;
    for (uint8_t i=1; i<=p; i++){
	z_poly[i] = err_poly[i];
	for (int j=i-1; j>=0; j--){
	    z_poly[i] = z_poly[i] ^ gf2m_mult(err_poly[j], syndroms[i-j-1]);
	}
	if (z_poly[i]){
	    z_poly_idx[i] = p_rs->index_of_alphas[ z_poly[i] ];
	}
    }

    /* find the error magnitude */
    for (int i=0; i<p; i++){
	uint8_t beta_inv_idx = n - beta_idx[i];
	uint8_t z = 1;
	uint8_t denom = 0;
	for (int j=1; j<=p; j++){
	    if (z_poly[j]){
		uint8_t idx = (z_poly_idx[j] + j*beta_inv_idx) % order;
		z = z ^ p_rs->alphas[ idx ];
	    }
	}

	for (int j=0; j<p; j++){
	    if (j == i){
		continue;
	    }
	    uint8_t d = 1 ^ p_rs->alphas [ EXPONET_ADD(beta_inv_idx,  beta_idx[j], order) ];
	    denom = EXPONET_ADD(denom, p_rs->index_of_alphas[d], order);
	}

	err_mag[i] = p_rs->alphas[EXPONET_ADD (p_rs->index_of_alphas[z], order-denom, order) ];
    }

    /* reuse syndroms scratch */
    memset(syndroms, 0, n);
    for (int i=0; i<p; i++){
	syndroms[ beta_idx[i] ] = err_mag[i];
    }

    for (int i=n-k; i<n; i++){
	out[i-(n-k)] =  in[i] ^ syndroms[i];
    }
    return k;
}


unsigned rs_decode_BM_with_erasure(rs_code         *p_rs
				  ,unsigned char   *in
				  ,int              in_sz
				  ,unsigned char   *erasures
				  ,int              nb_era
				  ,unsigned char   *out
				  ,int              out_sz
				  )
{
    const unsigned t = (p_rs->n - p_rs->k - nb_era) / 2;
    const unsigned n  = p_rs->n;
    const unsigned k  = p_rs->k;
    const uint8_t order = p_rs->mul_order;
    
    unsigned offset = 0;
    uint8_t	*syndroms      = p_rs->scratch + offset;  offset += n;   /* this will be reused, so n bytes */
    uint8_t	*erasure_poly  = p_rs->scratch + offset;  offset += nb_era+1;
    uint8_t	*era_poly_idx  = p_rs->scratch + offset;  offset += nb_era+1;    
    uint8_t	*mod_syndroms  = p_rs->scratch + offset;  offset += 2*t;
    uint8_t	*err_poly      = p_rs->scratch + offset;  offset += n-k+1;
    uint8_t	*err_coefs_idx = p_rs->scratch + offset;  offset += n-k+1;
    uint8_t	*beta_idx      = p_rs->scratch + offset;  offset += n-k+1;
    uint8_t	*z_poly	       = p_rs->scratch + offset;  offset += n-k+1;
    uint8_t	*z_poly_idx    = p_rs->scratch + offset;  offset += n-k+1;
    uint8_t	*err_mag       = p_rs->scratch + offset;  offset += n-k+1;
    uint8_t p = 0;

    
    assert(offset < 8*n);
    memset(erasure_poly, 0, nb_era+1);
    erasure_poly[0] = 1;
    era_poly_idx[0] = 0;
    
    for (int i=0; i<nb_era; i++){
	uint8_t expo = order - erasures[i];
	in[ erasures[i] ] = 0;
	for (int j=i+1; j>0; j--){
	    if (erasure_poly[j-1]){
		uint8_t e = EXPONET_ADD(era_poly_idx[j-1], expo, order);
		erasure_poly[j] = erasure_poly[j] ^ p_rs->alphas[e];
		era_poly_idx[j] = p_rs->index_of_alphas[ erasure_poly[j] ];
	    }
	}
    }
    
    rs_syndromes_calc(p_rs, in, syndroms, p_rs->scratch + offset);
    /* calculate the modified syndroms */
    uint8_t *syn_idx = z_poly;
    for (int i=0; i<(n-k); i++){
	syn_idx[i] = p_rs->index_of_alphas[ syndroms[i] ];
    }


    for (int i=0; i<n-k-nb_era; i++){
	mod_syndroms[i] = syndroms[i];
	for (int j=1; j<=nb_era; j++){
	    if (syndroms[i+j]){
		mod_syndroms[i] = mod_syndroms[i] ^ p_rs->alphas[ EXPONET_ADD (syn_idx[i+j], era_poly_idx[j], order) ];
	    }
	}
    }


    BM_iteration(p_rs, mod_syndroms, n-k-nb_era, err_poly, p_rs->scratch + offset);

    for (uint8_t i=0; i<=t; i++){
	err_coefs_idx[i] = p_rs->index_of_alphas[ err_poly[i] ];
    }
    
    for (int i=n-1; i>=0; i--){
	/* 1 test if i bit has errors, this is the same as testing
	 * if alpha^i is a root of err locator polynomial;
	 */
	uint8_t check = 1;
	for (uint8_t j=1; j<=t; j++){
	    if (err_poly[j]){
		err_coefs_idx[j] = EXPONET_ADD(err_coefs_idx[j], j, order);
		check = check ^ p_rs->alphas[ err_coefs_idx[j] ];
	    }
	}
	if (!check){
	    beta_idx[p]  = i;
	    p++;
	}
    }
    
    if (p > t){
	/* decoding failure */
	return 0;
    }
    else if (p==0 && nb_era == 0){
	memcpy(out, &in[n-k], k);
	return k;
    }

    /* at this point, we will add the erasures to the error locator */
    for (int i=p; i<p+nb_era; i++){
	uint8_t expo = erasures[i-p];
	err_poly[i+1] = 0;
	for (int j=i+1;j>0; j--){
	    if (err_poly[j-1]){
		uint8_t e = EXPONET_ADD(err_coefs_idx[j-1], expo, order);
		err_poly[j] = err_poly[j] ^ p_rs->alphas[e];
		err_coefs_idx[j] = p_rs->index_of_alphas[ err_poly[j] ];
	    }
	}
	beta_idx[i] = erasures[i-p];
    }
    p += nb_era;

    /* compute Z(x)  */
    z_poly[0] = 1;
    for (uint8_t i=1; i<=p; i++){
	z_poly[i] = err_poly[i];
	for (int j=i-1; j>=0; j--){
	    z_poly[i] = z_poly[i] ^ gf2m_mult(err_poly[j], syndroms[i-j-1]);
	}
	if (z_poly[i]){
	    z_poly_idx[i] = p_rs->index_of_alphas[ z_poly[i] ];
	}
    }

    /* find the error magnitude */
    for (int i=0; i<p; i++){
	uint8_t beta_inv_idx = n - beta_idx[i];
	uint8_t z = 1;
	uint8_t denom = 0;
	for (int j=1; j<=p; j++){
	    if (z_poly[j]){
		uint8_t idx = (z_poly_idx[j] + j*beta_inv_idx) % order;
		z = z ^ p_rs->alphas[ idx ];
	    }
	}
	if (z) {
	    for (int j=0; j<p; j++){
		if (j == i){
		    continue;
		}
		uint8_t d = 1 ^ p_rs->alphas [ EXPONET_ADD(beta_inv_idx,  beta_idx[j], order) ];
		denom = EXPONET_ADD(denom, p_rs->index_of_alphas[d], order);
	    }
	    
	    err_mag[i] = p_rs->alphas[EXPONET_ADD (p_rs->index_of_alphas[z], order-denom, order) ];
	}else{
	    err_mag[i] = 0;
	}
    }

    /* reuse syndroms scratch */
    memset(syndroms, 0, n);
    for (int i=0; i<p; i++){
	syndroms[ beta_idx[i] ] = err_mag[i];
    }

    for (int i=n-k; i<n; i++){
	out[i-(n-k)] =  in[i] ^ syndroms[i];
    }
    return k;
}
