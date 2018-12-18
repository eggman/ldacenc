/*
 * display ldac info
 *
 * url: https://github.com/eggman/ldacenc
 * licence: Apache License, Version 2.0
 */

#include <stdio.h>
#include "table.h"

AC g_ac0, g_ac1;
AB g_ab;

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

void dump_frame_header_ldac(unsigned char *pdata, int *p_loc)
{
    int syncword;

    syncword = read_bits(pdata, *p_loc + 0, 8);
    if (syncword != 0xAA) {
        printf("not found syncword\n");
        return;
    }

    printf("HEADER\n");
    printf("  SYNCWORD   %02X\n", syncword);
    printf("  SAMPLERATE %02X\n", read_bits(pdata, *p_loc +  8, 3));
    printf("  CHCONFIG   %02X\n", read_bits(pdata, *p_loc + 11, 2));
    printf("  FRAMELEN   %02X\n", read_bits(pdata, *p_loc + 13, 9));
    printf("  FRAMESTAT  %02X\n", read_bits(pdata, *p_loc + 22, 2));
    printf("\n");

    *p_loc += 24;
}

void dump_band_info_ldac(AB *p_ab, unsigned char *pdata, int *p_loc)
{
    printf("BANDINFO\n");
    printf("  NBAND      %02X\n", read_bits(pdata, *p_loc + 0, 4));
    p_ab->nbands = 2 + read_bits(pdata, *p_loc + 0, 4);
    p_ab->nqus = ga_nqus_ldac[p_ab->nbands];
    printf("  FLAG       %02X\n", read_bits(pdata, *p_loc + 4, 1));
    printf("\n");

    *p_loc += 5;
}

void dump_gradient_ldac(AB *p_ab, unsigned char *pdata, int *p_loc)
{
    p_ab->grad_mode = read_bits(pdata, *p_loc + 0, 2);

    printf("GRADIENT\n");
    printf("  GRADMODE   %02X\n", read_bits(pdata, *p_loc + 0, 2));
    if (p_ab->grad_mode == 0) {
        printf("  GRADQU0_L  %02X\n", read_bits(pdata, *p_loc +  2, 6));
        printf("  GRADQU0_H  %02X\n", read_bits(pdata, *p_loc +  8, 6));
        printf("  GRADOS_L   %02X\n", read_bits(pdata, *p_loc + 14, 5));
        printf("  GRADOS_H   %02X\n", read_bits(pdata, *p_loc + 19, 5));
        printf("  NADJQU     %02X\n", read_bits(pdata, *p_loc + 24, 5));
        p_ab->nadjqus = read_bits(pdata, *p_loc + 24, 5);

        int hqu = p_ab->nqus;
        int iqu;
        int grad_qu_l = read_bits(pdata, *p_loc +  2, 6);
        int grad_qu_h = read_bits(pdata, *p_loc +  8, 6) + 1;
        int grad_os_l = read_bits(pdata, *p_loc + 14, 5);
        int grad_os_h = read_bits(pdata, *p_loc + 19, 5);
        int tmp = grad_qu_h - grad_qu_l;
        int *p_grad = p_ab->a_grad;
        const unsigned char *p_t;

        /* Calculate Gradient Curve */
        for (iqu = 0; iqu < grad_qu_h; iqu++) {
            p_grad[iqu] = -grad_os_l;
        }
        for (iqu = grad_qu_h; iqu < hqu; iqu++) {
            p_grad[iqu] = -grad_os_h;
        }

        if (tmp > 0) {
            p_t = gaa_resamp_grad_ldac[tmp-1];

            tmp = grad_os_h - grad_os_l;
            if (tmp > 0) {
                tmp = tmp-1;
                for (iqu = grad_qu_l; iqu < grad_qu_h; iqu++) {
                    p_grad[iqu] -= ((*p_t++ * tmp) >> 8) + 1;
                }
            }
            else if (tmp < 0) {
                tmp = -tmp-1;
                for (iqu = grad_qu_l; iqu < grad_qu_h; iqu++) {
                    p_grad[iqu] += ((*p_t++ * tmp) >> 8) + 1;
                }
            }
        }
        printf("  a_grad   ");
        for (iqu = 0; iqu < hqu; iqu++ ) {
            printf(" %d", p_grad[iqu]);
        }
        printf("\n");
        *p_loc += 29;

    } else {
        printf("  GRADQU1    %02X\n", read_bits(pdata, *p_loc +  2, 5));
        printf("  GRADOS     %02X\n", read_bits(pdata, *p_loc +  7, 5));
        printf("  NADJQU     %02X\n", read_bits(pdata, *p_loc + 12, 5));
        p_ab->nadjqus = read_bits(pdata, *p_loc + 12, 5);
        *p_loc += 17;
    }
    printf("\n");
}

int decode_huffman(const CODES c, unsigned char *pdata, int start_pos)
{
    int pos = start_pos;
    int tmp = 0;
    int ret;

    /* pre-read  */
    if (c.codes_min_bits > 1) {
        tmp = read_bits(pdata, pos, c.codes_min_bits - 1);
        pos += c.codes_min_bits - 1;
    }

    /* check codes */
    for (int nbits = c.codes_min_bits; nbits <= c.codes_max_bits; nbits++) {
        tmp = tmp << 1;
        tmp += read_bit(pdata, pos);
        pos++;
        for (int i = 0; i < c.ncodes; i++) {
            if (c.p_tbl[i].len == nbits && c.p_tbl[i].code == tmp) {
                ret = i;
                goto found;
            }
        }
    }
    puts("huffman code not found");
    return -1;

found:
    return ret;
}

int dump_ldac_sfhuffman(AC *p_ac, unsigned char *pdata, int pos, const CODES c)
{
    int p, count = 1, idx;
    int dif[LDAC_MAXNQUS];

    printf("  DIF       ");

    if (p_ac->sfc_mode == 1 && p_ac->ich != 0) {
        count = 0;
    } else {
        count = 1;
    }

    for (p = pos; count < p_ac->p_ab->nqus;) {
        idx = decode_huffman(c, pdata, p);
        if (idx >= 0) {
            if (idx & (c.ncodes>>1)) {
                dif[count] = -c.ncodes + idx;
            } else {
                dif[count] = idx;
            }
            printf(" %2d", dif[count]);
            p += c.p_tbl[idx].len;
            count++;
        } else {
            break;
        }
    }
    printf("\n");

    /* compute idsf */
    if(p_ac->sfc_mode == 0) {
        int val0 = p_ac->a_idsf[0];
        int val1;
        const unsigned char *p_tbl;

        p_tbl = gaa_sfcwgt_ldac[p_ac->sfc_weight];
        p_ac->a_idsf[0] = val0 - p_tbl[0];

        for (int iqu = 1; iqu < p_ac->p_ab->nqus; iqu++) {
            val1 = dif[iqu] + val0;
            p_ac->a_idsf[iqu] = val1 - p_tbl[iqu];
            val0 = val1;
        }
    } else if (p_ac->sfc_mode == 1 && p_ac->ich != 0) {
        for (int iqu = 0; iqu < p_ac->p_ab->nqus; iqu++) {
            p_ac->a_idsf[iqu] = p_ac->p_ab->ap_ac[0]->a_idsf[iqu] + dif[iqu];
        }
    }

    printf("  a_idsf    ");
    for (int i = 0; i < p_ac->p_ab->nqus; i++) {
        printf(" %02X", p_ac->a_idsf[i]);
    }
    printf("\n");

    return p - pos;
}

void dump_scale_factor_0_ldac(AC *p_ac, unsigned char *pdata, int *p_loc)
{
    int sfc_bitlen, sfc_offset, hc_len;

    printf("  SFCBLEN    %02X\n", read_bits(pdata, *p_loc +  0, 2));
    sfc_bitlen = 3 + read_bits(pdata, *p_loc +  0, 2);
    p_ac->sfc_bitlen = sfc_bitlen;
    printf("  IDSF       %02X\n", read_bits(pdata, *p_loc +  2, 5));
    sfc_offset = read_bits(pdata, *p_loc +  2, 5);
    p_ac->sfc_offset = sfc_offset;
    printf("  SFCWTBL    %02X\n", read_bits(pdata, *p_loc +  7, 3));
    p_ac->sfc_weight = read_bits(pdata, *p_loc +  7, 3);

    printf("  VAL0       %02X\n", read_bits(pdata, *p_loc + 10, sfc_bitlen));
    p_ac->a_idsf[0] = read_bits(pdata, *p_loc + 10, sfc_bitlen) + sfc_offset;

    hc_len = dump_ldac_sfhuffman(p_ac, pdata, *p_loc + 10 + sfc_bitlen, codes0);
    *p_loc += 10 + sfc_bitlen + hc_len;
}

void dump_scale_factor_2_ldac(AC *p_ac, unsigned char *pdata, int *p_loc)
{
    int sfc_bitlen, hc_len;

    printf("  SFCBLEN    %02X\n", read_bits(pdata, *p_loc +  0, 2));
    sfc_bitlen = 2 + read_bits(pdata, *p_loc +  0, 2);

    hc_len = dump_ldac_sfhuffman(p_ac, pdata, *p_loc + 2, codes1);
    *p_loc += 2 + hc_len;
}

void dump_scale_factor_ldac(AC *p_ac, unsigned char *pdata, int *p_loc)
{
    int sfc_mode;

    printf("SCALEFACTOR\n");
    sfc_mode = read_bits(pdata, *p_loc +  0, 1);
    printf("  SFCMODE    %02X\n", sfc_mode);
    p_ac->sfc_mode = sfc_mode;
    *p_loc += 1;

    if (p_ac->ich == 0) {
        if (sfc_mode == 0) {
            dump_scale_factor_0_ldac(p_ac, pdata, p_loc);
        } else {
        }
    } else {
        if (sfc_mode == 0) {
            dump_scale_factor_0_ldac(p_ac, pdata, p_loc);
        } else {
            dump_scale_factor_2_ldac(p_ac, pdata, p_loc);
        }
    }
    printf("\n");
}

void dump_spectrum_ldac(AC *p_ac, unsigned char *pdata, int *p_loc)
{
    int i, iqu, idwl1, idwl2;
    int hqu = p_ac->p_ab->nqus;
    int nqus = hqu;
    int isp;
    int lsp, hsp;
    int nsps, wl;
    int *p_grad, *p_idsf, *p_idwl1, *p_idwl2, *p_tmp;

    printf("SPECTRUM\n");

    p_grad = p_ac->p_ab->a_grad;
    p_idsf = p_ac->a_idsf;
    p_idwl1 = p_ac->a_idwl1;
    p_idwl2 = p_ac->a_idwl2;

    /* grad_mode == 0 */
    for (iqu = 0; iqu < hqu; iqu++) {
        for (iqu = 0; iqu < hqu; iqu++) {
            idwl1 = p_idsf[iqu] + p_grad[iqu];
            if (idwl1 < LDAC_MINIDWL1) {
                idwl1 = LDAC_MINIDWL1;
            }
            idwl2 = 0;
            if (idwl1 > LDAC_MAXIDWL1) {
                idwl2 = idwl1 - LDAC_MAXIDWL1;
                if (idwl2 > LDAC_MAXIDWL2) {
                    idwl2 = LDAC_MAXIDWL2;
                }
                idwl1 = LDAC_MAXIDWL1;
            }
            p_idwl1[iqu] = idwl1;
            p_idwl2[iqu] = idwl2;
        }
    }

    printf("  a_idwl1 a ");
    for (iqu = 0; iqu < nqus; iqu++) {
        printf(" %02X", p_idwl1[iqu]);
    }
    printf("\n");

    /* adjust */
    int nadjqus = p_ac->p_ab->nadjqus;
	p_tmp = p_ac->a_idwl1;
    for (iqu = 0; iqu < nqus; iqu++) {
        idwl1 = p_tmp[iqu];
        if (iqu < nadjqus) {
            idwl1++;
        }
        idwl2 = 0;
        if (idwl1 > LDAC_MAXIDWL1) {
            idwl2 = idwl1 - LDAC_MAXIDWL1;
            if (idwl2 > LDAC_MAXIDWL2) {
                idwl2 = LDAC_MAXIDWL2;
            }
            idwl1 = LDAC_MAXIDWL1;
        }
        p_idwl1[iqu] = idwl1;
        p_idwl2[iqu] = idwl2;
    }

    printf("  a_idwl1 b ");
    for (iqu = 0; iqu < nqus; iqu++) {
        printf(" %02X", p_idwl1[iqu]);
    }
    printf("\n");

    /* wl */
    int p_wl[LDAC_MAXLSU];
    int j = 0;
    printf("  wl        ");
    for (iqu = 0; iqu < nqus; iqu++) {
        lsp = ga_isp_ldac[iqu];
        hsp = ga_isp_ldac[iqu+1];
        nsps = ga_nsps_ldac[iqu];
        idwl1 = p_ac->a_idwl1[iqu];
        wl = ga_wl_ldac[idwl1];

        if (idwl1 == 1) {
            isp = lsp;

            if (nsps == 2) {
                printf(" %2d", LDAC_2DIMSPECBITS);
                p_wl[j] = LDAC_2DIMSPECBITS;
                j++;
            } else {
                for (i = 0; i < nsps>>2; i++, isp+=4) {
                   printf(" %2d", LDAC_4DIMSPECBITS);
                   p_wl[j] = LDAC_4DIMSPECBITS;
                   j++;
                }
            }
        } else {
            for (isp = lsp; isp < hsp; isp++) {
                printf(" %2d", wl);
                p_wl[j] = wl;
                j++;
            }
        }
    }
    printf("\n");

    // spectrum
    int nsp = j;

    printf(" coded spec");
    for (i = 0; i < nsp; i++) {
        printf(" %03X", read_bits(pdata, *p_loc, p_wl[i]));
        *p_loc += p_wl[i];
    }
    printf("\n\n");
}

void dump_residual_ldac(AC *p_ac, unsigned char *pdata, int *p_loc)
{
    int iqu, isp;
    int lsp, hsp;
    int nqus = p_ac->p_ab->nqus;
    int idwl2, wl;
    int *p_idwl2;

    p_idwl2 = p_ac->a_idwl2;

    printf("RESIDUAL\n");

    printf("  a_idwl2   ");
    for (iqu = 0; iqu < nqus; iqu++) {
        printf(" %02X", p_idwl2[iqu]);
    }
    printf("\n");

    printf(" coded spec");
    for (iqu = 0; iqu < nqus; iqu++) {
        idwl2 = p_idwl2[iqu];

        if (idwl2 > 0) {
            lsp = ga_isp_ldac[iqu];
            hsp = ga_isp_ldac[iqu+1];
            wl = ga_wl_ldac[idwl2];

            for (isp = lsp; isp < hsp; isp++) {
                printf(" %03X", read_bits(pdata, *p_loc, wl));
                *p_loc += wl;
            }
        }
    }
    printf("\n\n");
}

void dump_byte_alignment_ldac(unsigned char *pdata, int *p_loc)
{
    int nbits_padding;

    nbits_padding = ((*p_loc + LDAC_BYTESIZE - 1) / LDAC_BYTESIZE) * LDAC_BYTESIZE - *p_loc;

    if (nbits_padding > 0) {
        printf("PADDING %d bits\n\n", nbits_padding);
        *p_loc += nbits_padding;
    }

    return;
}

int main(int argc, char *argv[])
{
    int pos, *p_loc;
    unsigned char ldac[1024];
    FILE *infp;
    AB *p_ab;
    AC *p_ac;

    /* initialize little */
    p_ab = &g_ab;
    p_ab->ap_ac[0] = &g_ac0;
    p_ab->ap_ac[1] = &g_ac1;
    p_ab->ap_ac[0]->ich = 0;
    p_ab->ap_ac[0]->p_ab = p_ab;
    p_ab->ap_ac[1]->ich = 1;
    p_ab->ap_ac[1]->p_ab = p_ab;

    pos = 0;
    p_loc = &pos;

    if ((infp = fopen(argv[1], "r"))==NULL) {
        printf("can't open input file\n");
        return -1;
    }

    fread(ldac, 660, 1, infp);

    dump_frame_header_ldac(ldac, p_loc);
    dump_band_info_ldac(p_ab, ldac, p_loc);
    dump_gradient_ldac(p_ab, ldac, p_loc);

    p_ac = p_ab->ap_ac[0];
    dump_scale_factor_ldac(p_ac, ldac, p_loc);
    dump_spectrum_ldac(p_ac, ldac, p_loc);
    dump_residual_ldac(p_ac, ldac, p_loc);
    p_ac = p_ab->ap_ac[1];
    dump_scale_factor_ldac(p_ac, ldac, p_loc);
    dump_spectrum_ldac(p_ac, ldac, p_loc);
    dump_residual_ldac(p_ac, ldac, p_loc);

    dump_byte_alignment_ldac(ldac, p_loc);

    dump_frame_header_ldac(ldac, p_loc);

    return 0;
}

