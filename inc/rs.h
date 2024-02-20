#pragma once
#include <stdio.h>
#include <stdlib.h>

typedef struct rs255_encoder_s  rs255_encoder;


rs255_encoder* rs255_encoder_construct(int k);
void rs255_encoder_destroy(rs255_encoder *p);
unsigned rs255_encode(rs255_encoder *p_enc
                     ,unsigned char      *in
                     ,int              in_sz
                     ,unsigned char      *out
                     ,int              out_sz
                     );
