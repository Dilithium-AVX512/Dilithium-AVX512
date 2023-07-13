#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "../align.h"
#include "../params.h"
#include "../sign.h"
#include "../packing.h"
#include "../polyvec.h"
#include "../poly.h"
#include "../randombytes.h"
#include "../symmetric.h"
#include "../fips202.h"
#ifdef DILITHIUM_USE_AES
#include "../aes256ctr.h"
#endif

#define NTESTS 100
#define MLEN 59

int main()
{
    size_t i = 0, j = 0;
    int ret;
    size_t mlen, smlen;
    ;
    uint8_t m[MLEN + CRYPTO_BYTES];
    uint8_t sm[MLEN + CRYPTO_BYTES];
    uint8_t sm_f[MLEN + CRYPTO_BYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];

    FILE *fp = fopen("sign_kat.txt", "r+");
    char prefix[5];
    char prefi[4];
    char pref[5];

    while (j < NTESTS)
    {
        printf("%d: ", j);
        fgets(prefix, sizeof(prefix), fp);
        i = 0;
        // printf("sk = ");
        while (fscanf(fp, "%u", &sk[i]) && i < CRYPTO_SECRETKEYBYTES)
        {

            // printf("%u ", sk[i]);
            // if (i % 8 == 7)
            //     printf("\n");
            i++;
        }
        // printf("\n");

        fgets(prefi, sizeof(prefi), fp);
        // printf("m = ");
        i = 0;
        while (fscanf(fp, "%u", &m[i]) && i < MLEN)
        {
            // printf("%u ", m[i]);
            // if (i % 8 == 7)
            //     printf("\n");
            i++;
        }
        // printf("\n");
        // printf("\n");

        crypto_sign(sm, &smlen, m, MLEN, sk);

        fgets(pref, sizeof(pref), fp);
        // printf("sm_f = ");
        i = 0;
        while (fscanf(fp, "%u", &sm_f[i]) && i < smlen)
        {
            // printf("%u ", sm_f[i]);
            // if (i % 8 == 7)
            //     printf("\n");
            i++;
        }
        // printf("\n");

        i = 0;
        if ((memcmp(sm, sm_f, smlen)) != 0 )
        {

            printf("fail\n");
            return 0;
        }
        else
            printf("ok\n");

        j++;
    }

    return 0;
}
