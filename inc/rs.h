#pragma once
#include <stdio.h>
#include <stdlib.h>

typedef struct rs_code_s  rs_code;

rs_code* rs_code_construct(int n, int k);
void rs_code_destroy(rs_code *p);
unsigned rs_encode(rs_code *p_enc
		  ,unsigned char      *in
		  ,int              in_sz
		  ,unsigned char      *out
		  ,int              out_sz
		  );

int rs_get_genpoly(rs_code *enc, uint8_t* out, int out_sz);
unsigned rs_encode_sys(rs_code *p_enc
		      ,unsigned char      *in
		      ,int              in_sz
		      ,unsigned char      *out
		      ,int              out_sz
		      );

unsigned rs_decode_ge(rs_code         *p_rs
		     ,unsigned char   *in
		     ,int              in_sz
		     ,unsigned char   *out
		     ,int              out_sz
		     );


unsigned rs_decode_BM(rs_code         *p_rs
		     ,unsigned char   *in
		     ,int              in_sz
		     ,unsigned char   *out
		     ,int              out_sz
		     );
