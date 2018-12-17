/*
 * display ldac info
 *
 * url: https://github.com/eggman/ldacenc
 * licence: Apache License, Version 2.0
 */

#include <stdio.h>

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

int dump_ldac_scalefactor(unsigned char *pdata, int pos)
{
    int sfc_bitlen;
    printf("SCALEFACTOR\n");
    printf("  SFCMODE    %02X\n", read_bits(pdata, pos +  0, 1));
    printf("  SFCBLEN    %02X\n", read_bits(pdata, pos +  1, 2));
    sfc_bitlen = 3 + read_bits(pdata, pos +  1, 2);
    printf("  IDSF       %02X\n", read_bits(pdata, pos +  3, 5));
    printf("  SFCWTBL    %02X\n", read_bits(pdata, pos +  8, 3));

    printf("  VAL0       %02X\n", read_bits(pdata, pos + 11, sfc_bitlen));

    // todo decode huffman

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

