/*polynomial operations under 2^m finite field, where m <= 8 */
#include "gf2_8.h"


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
