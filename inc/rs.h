#pragma once
#include <stdio.h>
#include <stdlib.h>

typedef struct rs_encoder_s  rs_encoder;

rs_encoder* rs_encoder_construct(int n, int k);
void rs_encoder_destroy(rs_encoder *p);
unsigned rs_encode(rs_encoder *p_enc
		  ,unsigned char      *in
		  ,int              in_sz
		  ,unsigned char      *out
		  ,int              out_sz
		  );

int rs_get_genpoly(rs_encoder *enc, uint8_t* out, int out_sz);
