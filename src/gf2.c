#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "gf2.h"

/* convolve the coefficient of f and g, deg(f) and deg(g) < 8 */
uint16_t gf2_poly_mult(uint8_t f, uint8_t g)
{
    uint8_t len_g = deg_plus_one(g);
    uint16_t r = 0;
    uint16_t x = f;

    for (uint8_t i=0; i<len_g; i++){
        uint16_t mask = (g & 0x01) ? 0xFFFF : 0x0000;
        r ^= (x & mask);
        x = x << 1;
        g = g >> 1;
    }

    return r;
}


/* find prime polynomials upto a byte, deg(f) <= 8, we eleminiate all
   polynoimials can be produced by g(x)*h(x).
*/
unsigned find_prime_poly(int deg, unsigned* out, int out_sz)
{
    uint8_t *cand;
    int len = 0;
    if (deg > 8){
        fprintf(stderr, "not implemented yet, deg has to be <= 8");
        return 0;
    }

    cand = calloc(512, 1);

    for (unsigned h=1; h<deg; h++){
        unsigned g=deg - h;
        if (g < h){
            break;
        }
        //enumerate all the polys
        for (unsigned i=(1<<h); i < (1<<(h+1)); i++){
            for (unsigned j=(1<<g); j < (1<<(g+1)); j++){
                uint16_t r = gf2_poly_mult(i, j);

                assert(r < 512);
                cand[r] = 1;
            }
        }
    }

    for(unsigned i= 1<<deg; i< 1<<(deg+1); i++){
        if(!cand[i]){
            out[len++] = i;
        }
    }

    free(cand);
    return len;
}

/* returns [ a(x) * b(x) ] mod g(x) where g(x) is a prime polynomial
 * all degrees of a(x) b(x) and g(x) is less then 8
 */
uint8_t gf2_8_mult_mod(uint8_t a, uint8_t b, unsigned g)
{
    uint16_t ab = gf2_poly_mult(a, b);
    uint8_t d_ab = deg_plus_one(ab);
    uint8_t d_g = deg_plus_one(g);

    while(d_ab >= d_g){
        uint8_t q = d_ab - d_g;
        ab = ab ^ (g << q);
        d_ab = deg_plus_one(ab);
    }
    return ab & 0xFF;
}
        

/* find the primitive element given the field elements are in F_g[x]
 *
 * we will brutal search for each nonzero elements to see if it becomes
 * one for any power before q-1, and will return the first one that is
 * primitive
 */
uint8_t find_primitive_element(unsigned g)
{
    unsigned dg = deg_plus_one(g) -1;
    int b_found = 0;
    assert(dg <= 8);

    for (unsigned i=0;i<(1<<dg);i++){
        uint8_t a = i & 0xFF;
        uint8_t r = a;
        unsigned j = 2;
        for (j=2; j<(1<<dg)-1; j++){
            r = gf2_8_mult_mod(r, a, g);
            //printf("%d^%d = %d \n ", a, j, r);
            if (r == a){
                //printf(" %d, %d is not a primitive\n", j, a);
                break;
            }
        }
        if (j == (1<<dg) -1){
            r =  gf2_8_mult_mod(r, a, g);
            assert(r == 1);
            return a;
        }
    }
    return 0;
}
        
    

    

    
