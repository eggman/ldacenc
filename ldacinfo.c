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

    return 24;
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

    return 0;
}

