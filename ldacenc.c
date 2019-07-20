/*
 * ldac enocoder
 *
 * url: https://github.com/eggman/ldacenc
 * licence: Apache License, Version 2.0
 */

#include <stdio.h>
#include "ldacBT.h"

int main(int argc, char *argv[])
{
    HANDLE_LDAC_BT h;
    int ret, err, pcm_used, ldac_stream_size, ldac_frame_num;
    unsigned char pcm[128*2*2];
    unsigned char ldac[1024];
    FILE *infp, *outfp;
    unsigned char header[44];

    if ((infp = fopen(argv[1], "r"))==NULL) {
        printf("can't open input file\n");
        return -1;
    }

    if ((outfp = fopen(argv[2], "w"))==NULL) {
        printf("can't open output file\n");
        return -1;
    }

    //skip wav header
    if (1 != fread(header, 44, 1, infp) ) {
        return -1;
    }
    
    h = ldacBT_get_handle();
    ret = ldacBT_init_handle_encode(h, 679, LDACBT_EQMID_MQ, LDACBT_CHANNEL_MODE_STEREO, LDACBT_SMPL_FMT_S16, 48000);

    while(fread(pcm, 512, 1, infp))
    {
        ret = ldacBT_encode(h, &pcm[0], &pcm_used, ldac, &ldac_stream_size, &ldac_frame_num);
        printf("  pcm_used=%d ldac_stream_size=%d ldac_frame_num=%d\n", pcm_used, ldac_stream_size, ldac_frame_num);
        fwrite(ldac, 1, ldac_stream_size, outfp);
    }

    ldacBT_close_handle(h);

    return 0;
}

