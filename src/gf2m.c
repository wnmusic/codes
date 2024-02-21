/*polynomial operations under 2^m finite field, where m <= 8 */
#include "gf2m.h"

unsigned num_ele_mutiplicity;
const uint8_t *alphas = NULL;
const uint8_t *index_of_alphas = NULL;

/* multiply f(x) with g(x), deg[f(x)] < df, deg[g(x)]  < dg
 * the coeffients are from Fq, where q=2^m
 * by convention, f array stores [f0, f1, .... f_{df-1}] 
 */
int gf2m_poly_mult(uint8_t* f, int df, uint8_t* g, int dg, uint8_t *out, int out_sz)
{
    int dout = df + dg - 1;
    if (out_sz < dout){
        return 0;
    }
    
    for (int i=0; i<dout; i++){
        out[i] = 0;
        int start = i-dg +1 > 0 ? i-dg+1 : 0;
        for (int j=start; j<=i && j<df ;j++){
            out[i] = gf2m_add(out[i], gf2m_mult(f[j], g[i-j]));
        }
    }
    return dout;
}

/* divide a polynomial f(x) by g(x), the quotient will be stored in q(x) and the
 * reminder will be r(x), deg[r(x)] < deg[g(x)] and will be stored in f
 */
void gf2m_poly_div(uint8_t* fr, int dfr,
                  uint8_t* g, int dg,
                  uint8_t *out, int out_sz,
                  int *dq, int *dr)
{

    memset(out, 0, sizeof(uint8_t) * out_sz);
    *dq = dfr - dg + 1;
    while(dfr >= dg){
        uint8_t m = gf2m_div(fr[dfr-1], g[dg-1]);
        int new_df = 0;
        
        out[dfr - dg] = m;
        for (unsigned i=dfr -dg, j=0; i<dfr; i++,j++){
            fr[i] = gf2m_add(fr[i], gf2m_mult(m, g[j]));
        }

        for (unsigned i=0; i<dfr; i++){
            if (fr[i]){
                new_df = i+1;
            }
        }
        dfr = new_df;
    }
    *dr = dfr;
}


/* the matrix inverse on GF(2^m). matrix A will be square,
 * After operations, A will be the ideintify matrix,
 * and b contains the result 
 */
void gf2m_gaussian_elemination(uint8_t *A, int dim_A1, int dim_A2, uint8_t *b, int dim)
{
    assert(dim_A1 == dim_A2);
    assert(dim_A2 == dim);

    (void)dim_A1;
    (void)dim_A2;

    /* forward round */
    for (unsigned i=0; i<dim-1; i++){
        for (unsigned j=i+1; j<dim; j++){
            
            uint8_t m = gf2m_div(A[j*dim+i], A[i*dim+i]);
            
            for (unsigned k=i;k<dim; k++){
                A[j*dim+k] = gf2m_add(A[j*dim+k], gf2m_mult(m, A[i*dim + k]));
            }

            b[j] = gf2m_add(b[j], gf2m_mult(m, b[i]));
        }
    }

    /* bottom up */
    for (int i=dim-1; i> 0; i--){
        for (int j=i-1; j>=0; j--){
            uint8_t m = gf2m_div(A[j*dim+i], A[i*dim+i]);
            b[j] = gf2m_add(b[j], gf2m_mult(m, b[i]));
            A[j*dim + i] = 0;
        }
    }

    for (unsigned i=0; i<dim; i++){
        b[i] = gf2m_div(b[i], A[i*dim +i]);
        A[i*dim + i] = 1;
    }

}
    
