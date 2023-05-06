#pragma once

#ifdef _DEBUG
#include <stdio.h>
#include <string.h>
static FILE *fp;

typedef struct{
    char name[64];
    unsigned ele_sz; /* sizeof each element */
    unsigned nb_ele; /* how many elements per frame */
    char ctype[16]; /* float or int ... */
}trace_info;

#define MAX_NUM_TRACES   (64)
static FILE *fp = NULL;
static trace_info trace_header[MAX_NUM_TRACES] = {};
static unsigned gb_num_traces = 0;


void trace_open()
{
    char fname[128];
    int l;
    sprintf(fname, "%s", __FILE__);
    l = strlen(fname);
    fname[l-1] = 'd';
    fname[l] = 'a';
    fname[l+1] = 't';
    fname[l+2] = 0;
        
    fp = fopen("trace.dat", "wb");
    memset(trace_header, 0, sizeof(trace_header));
    
    fseek(fp, sizeof(trace_header), SEEK_SET);

}

void trace_close()
{
    fseek(fp, 0L, SEEK_SET);
    fwrite(&trace_header[0], sizeof(trace_header), 1, fp);
    fseek(fp, 0L, SEEK_END);
    fclose(fp);
}

#define TRY_ADD_TRACE_HEADER(_name, _ctype, _num) do {\
    unsigned i;                                      \
    for (i=0; i<gb_num_traces; i++){                 \
        if (strncmp(_name, trace_header[i].name, 64)== 0 ){  \
            break;                                   \
        }                                            \
    }                                                \
    if (i == gb_num_traces){                         \
        strncpy(trace_header[gb_num_traces].name, _name, 64);      \
        strncpy(trace_header[gb_num_traces].ctype, #_ctype, 16);   \
        trace_header[gb_num_traces].ele_sz = sizeof(_ctype);\
        trace_header[gb_num_traces].nb_ele = _num;     \
        gb_num_traces++;                             \
     }                                            \
}while(0)

#define DBG_TRACE_OPEN()         trace_open()
#define DBG_TRACE_CLOSE()        trace_close()

#define DBG_TRACE_VALUE(name, ctype, val)                 do{\
    TRY_ADD_TRACE_HEADER(name, ctype, 1);                               \
    fwrite(&val, sizeof(ctype), 1, fp);                                 \
    }while(0)

#define DBG_TRACE_ARRAY_1D(name, ctype, pval, dim)  do{\
    TRY_ADD_TRACE_HEADER(name, ctype, dim);                               \
    fwrite(pval, sizeof(ctype), dim, fp);                                 \
    }while(0)



#else
#define DBG_TRACE_OPEN();
#define DBG_TRACE_CLOSE();

#define DBG_TRACE_VALUE(name, ctype, val) ;
#define DBG_TRACE_ARRAY_1D(name, ctype, pval, dim);

#endif
