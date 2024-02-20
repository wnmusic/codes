#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

uint16_t gf2_poly_mult(uint8_t f, uint8_t g);
unsigned find_prime_poly(int deg, unsigned* out, int out_sz);
uint8_t gf2_8_mult_mod(uint8_t a, uint8_t b, unsigned g);
uint8_t find_primitive_element(unsigned g);

int gf2m_poly_mult(uint8_t* f, int df, uint8_t* g, int dg, uint8_t *out, int out_sz);
void gf2m_poly_div(uint8_t* f, int df,
                  uint8_t* g, int dg,
                  uint8_t *out, int out_sz,
                  int *dq, int *dr);

