/*
 * display ldac info
 *
 * url: https://github.com/eggman/ldacenc
 * licence: Apache License, Version 2.0
 */

#include <stdio.h>

#define DECLFUNC static
#define LDAC_NSFCWTBL          8
#define LDAC_MAXNQUS          34
DECLFUNC const unsigned char gaa_sfcwgt_ldac[LDAC_NSFCWTBL][LDAC_MAXNQUS] = {
{
     1,  0,  0,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,
     3,  3,  3,  3,  3,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,
},
{
     0,  1,  1,  2,  3,  4,  4,  4,  4,  5,  6,  6,  6,  6,  6,  7,
     7,  7,  7,  7,  7,  7,  8,  8,  8,  9, 10, 10, 11, 11, 12, 12, 12, 12,
},
{
     0,  1,  1,  2,  3,  3,  3,  3,  3,  4,  4,  5,  5,  5,  5,  5,
     5,  5,  5,  5,  5,  5,  6,  6,  6,  7,  8,  9,  9, 10, 10, 11, 11, 11,
},
{
     0,  1,  3,  4,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,
     7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  9,  9,  9, 10, 10, 10, 10,
},
{
     0,  1,  3,  4,  5,  5,  6,  7,  7,  8,  8,  9,  9, 10, 10, 10,
    10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13,
},
{
     1,  0,  1,  2,  2,  3,  3,  4,  4,  5,  6,  7,  7,  8,  8,  8,
     9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11,
},
{
     0,  0,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,
     4,  4,  4,  4,  4,  4,  4,  5,  5,  6,  7,  7,  7,  8,  9,  9,  9,  9,
},
{
     0,  0,  1,  2,  3,  4,  4,  5,  5,  6,  7,  7,  8,  8,  8,  8,
     9,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12,
},
};



typedef struct {
    unsigned char code;
    unsigned char len;
} CODEBOOK;

typedef struct {
    const CODEBOOK *p_tbl;
    unsigned char ncodes;
    unsigned char codes_min_bits;
    unsigned char codes_max_bits;
} CODES;

static const CODEBOOK codebook[8] = {
    { 0, 2}, { 1, 2}, {14, 4}, {62, 6},
    {63, 6}, {30, 5}, { 6, 3}, { 2, 2},
};

static const CODES codes = {
    codebook, 8, 2, 6
};

static int g_sfc_weight;
int g_a_idsf[LDAC_MAXNQUS];

int read_bit(unsigned char *pdata, int pos)
{
    int bytepos = pos / 8;
    int bitpos  = pos % 8;
    return (pdata[bytepos] & (1 << (7 - bitpos))) ? 1 : 0;
}

int read_bits(unsigned char *pdata, int pos, int nbits)
{
    int tmp = 0;
    int p = pos;

	for (int i = 1; i <= nbits; i++) {
        tmp = tmp << 1;
	    tmp += read_bit(pdata, p);
        p++;
    }
    return tmp;
}

int dump_ldac_header(unsigned char *pdata, int pos)
{
    int syncword;

    syncword = read_bits(pdata, pos + 0, 8);
    if (syncword != 0xAA) {
        printf("not found syncword\n");
        return -1;
    }

    printf("HEADER\n");
    printf("  SYNCWORD   %02X\n", syncword);
    printf("  SAMPLERATE %02X\n", read_bits(pdata, pos +  8, 3));
    printf("  CHCONFIG   %02X\n", read_bits(pdata, pos + 11, 2));
    printf("  FRAMELEN   %02X\n", read_bits(pdata, pos + 13, 9));
    printf("  FRAMESTAT  %02X\n", read_bits(pdata, pos + 22, 2));
    printf("\n");

    return pos + 24;
}

int dump_ldac_bandinfo(unsigned char *pdata, int pos)
{
    printf("BANDINFO\n");
    printf("  NBAND      %02X\n", read_bits(pdata, pos + 0, 4));
    printf("  FLAG       %02X\n", read_bits(pdata, pos + 4, 1));
    printf("\n");

    return pos + 5;
}

int dump_ldac_gradient(unsigned char *pdata, int pos)
{
    int new_pos = pos;
    int grad_mode = read_bits(pdata, pos + 0, 2);

    printf("GRADIENT\n");
        printf("  GRADMODE   %02X\n", read_bits(pdata, pos + 0, 2));
    if (grad_mode == 0) {
        printf("  GRADQU0    %02X\n", read_bits(pdata, pos +  2, 6));
        printf("  GRADQU0    %02X\n", read_bits(pdata, pos +  8, 6));
        printf("  GRADQS     %02X\n", read_bits(pdata, pos + 14, 5));
        printf("  GRADQS     %02X\n", read_bits(pdata, pos + 19, 5));
        printf("  NADJQU     %02X\n", read_bits(pdata, pos + 24, 5));
        new_pos = pos + 29;
    } else if (grad_mode == 1 || grad_mode == 2 || grad_mode == 3) {
        printf("  GRADQU1    %02X\n", read_bits(pdata, pos +  2, 5));
        printf("  GRADQS     %02X\n", read_bits(pdata, pos +  7, 5));
        printf("  NADJQU     %02X\n", read_bits(pdata, pos + 12, 5));
        new_pos = pos + 17;
    } else {
        printf("grad mode error\n");
        return -1;
    }
    printf("\n");

    return new_pos;
}

int decode_huffman(const CODES c, unsigned char *pdata, int start_pos)
{
    int pos = start_pos;
    int tmp;
    int i;

    /* read 1st bit */
    tmp =  read_bit(pdata, pos);

    /* check codes */
    for (int nbits = c.codes_min_bits; nbits <= c.codes_max_bits; nbits++) {
        pos++;
        tmp = tmp << 1;
        tmp += read_bit(pdata, pos);
        for (i=0; i < c.ncodes; i++) {
            if (c.p_tbl[i].code == tmp) {
                goto found;
            }
        }
    }
    puts("not found");
    return -1;
found:

    return i;
}

int dump_ldac_sfhuffman(unsigned char *pdata, int pos)
{
    int p, count, idx;
    int dif[LDAC_MAXNQUS];

    for (p = pos, count = 1; count < 24;) {
        idx = decode_huffman(codes, pdata, p);
        if (idx >= 0) {
            if (idx <= 4) {
                dif[count] = idx;
            } else {
                dif[count] = -8 + idx;
            }
            printf("  DIFF %2d    %2d (%d, %X)\n",
                   count, dif[count], codes.p_tbl[idx].len, codes.p_tbl[idx].code);
            p += codes.p_tbl[idx].len;
            count++;
        } else {
            break;
        }
    }

    //dump p_idsf
    int val0 = g_a_idsf[0];
    int val1;
    const unsigned char *p_tbl;

    p_tbl = gaa_sfcwgt_ldac[g_sfc_weight];
    g_a_idsf[0] = val0 - p_tbl[0];

    for (int iqu = 1; iqu < 24; iqu++) {
        val1 = dif[iqu] + val0;
        g_a_idsf[iqu] = val1 - p_tbl[iqu];
        val0 = val1;
    }

    for (int i = 0; i < 24 /* LDAC_MAXNQUS */; i++) {
        printf("  idsf[%2d]   %02X\n", i, g_a_idsf[i]);
    }

    return p - pos;
}

int dump_ldac_scalefactor(unsigned char *pdata, int pos)
{
    int sfc_bitlen, sfc_offset;
    printf("SCALEFACTOR\n");
    printf("  SFCMODE    %02X\n", read_bits(pdata, pos +  0, 1));
    printf("  SFCBLEN    %02X\n", read_bits(pdata, pos +  1, 2));
    sfc_bitlen = 3 + read_bits(pdata, pos +  1, 2);
    printf("  IDSF       %02X\n", read_bits(pdata, pos +  3, 5));
    sfc_offset = read_bits(pdata, pos +  3, 5);
    printf("  SFCWTBL    %02X\n", read_bits(pdata, pos +  8, 3));
    g_sfc_weight = read_bits(pdata, pos +  8, 3);

    printf("  VAL0       %02X\n", read_bits(pdata, pos + 11, sfc_bitlen));
    g_a_idsf[0] = read_bits(pdata, pos + 11, sfc_bitlen) + sfc_offset;

    dump_ldac_sfhuffman(pdata, pos + 11 + sfc_bitlen);

    return pos + 11 + sfc_bitlen;
}

int main(int argc, char *argv[])
{
    int pos;
    unsigned char ldac[1024];
    FILE *infp;

    if ((infp = fopen(argv[1], "r"))==NULL) {
        printf("can't open input file\n");
        return -1;
    }

    fread(ldac, 660, 1, infp);

    pos = dump_ldac_header(ldac, 0);
    pos = dump_ldac_bandinfo(ldac, pos);
    pos = dump_ldac_gradient(ldac, pos);
    pos = dump_ldac_scalefactor(ldac, pos);

    return 0;
}

