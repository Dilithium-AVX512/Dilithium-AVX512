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

int main()
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t pk_f[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk_f[CRYPTO_SECRETKEYBYTES];
    unsigned int i,j=0;
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
    const uint8_t *rho, *rhoprime, *key;
#ifdef DILITHIUM_USE_AES
    uint64_t nonce;
    aes256ctr_ctx aesctx;
    polyvecl rowbuf[1];
#else
    polyvecl rowbuf[2], mat[K];
#endif
    polyvecl s1, *row = rowbuf;
    polyveck s2;
    poly t1, t0;
#ifndef DILITHIUM_USE_AES
#if K == 6 && L == 5
    poly t2, t3, t4, t5;
#endif
#endif

    FILE *fp = fopen("keygen_kat.txt", "r+");
    char prefix[10];
    char prefi[5];
    i = 0;
    while (j < NTESTS)
    {
        printf("%d: ",j);
        fgets(prefix, sizeof(prefix), fp);
        i=0;
        // printf("seedbuf = ");
        while(fscanf(fp,"%u",&seedbuf[i]) && i<32){
            
            // printf("%u ",seedbuf[i]);
            // if(i%8==7) printf("\n");
            i++;
        }
        // printf("\n");

        fgets(prefix, sizeof(prefi), fp);
        // printf("pk = ");
        i=0;
        while(fscanf(fp,"%u",&pk_f[i]) && i<CRYPTO_PUBLICKEYBYTES){
            // printf("%u ",pk_f[i]);
            // if(i%8==7) printf("\n");
            i++;
        }
        // printf("\n");
        fgets(prefix, sizeof(prefi), fp);
        // printf("sk = \n");
        i=0;
        while(fscanf(fp,"%u",&sk_f[i]) && i<CRYPTO_SECRETKEYBYTES){
            // printf("%u ",sk_f[i]);
            // if(i%8==7) printf("\n");
            i++;
        }



        shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);
        rho = seedbuf;
        rhoprime = rho + SEEDBYTES;
        key = rhoprime + CRHBYTES;

        /* Store rho, key */
        memcpy(pk, rho, SEEDBYTES);
        memcpy(sk, rho, SEEDBYTES);
        memcpy(sk + SEEDBYTES, key, SEEDBYTES);

        /* Sample short vectors s1 and s2 */
#ifdef DILITHIUM_USE_AES
        aes256ctr_init(&aesctx, rhoprime, 0);
        for (i = 0; i < L; ++i)
        {
            nonce = i;
            aesctx.n = _mm_loadl_epi64((__m128i *)&nonce);
            poly_uniform_eta_preinit(&s1.vec[i], &aesctx);
        }
        for (i = 0; i < K; ++i)
        {
            nonce = L + i;
            aesctx.n = _mm_loadl_epi64((__m128i *)&nonce);
            poly_uniform_eta_preinit(&s2.vec[i], &aesctx);
        }
#elif K == 4 && L == 4
        poly_uniform_eta_8x(&s1.vec[0], &s1.vec[1], &s1.vec[2], &s1.vec[3], &s2.vec[0], &s2.vec[1], &s2.vec[2], &s2.vec[3], rhoprime, 0, 1, 2, 3, 4, 5, 6, 7);
#elif K == 6 && L == 5
        poly_uniform_eta_8x(&s1.vec[0], &s1.vec[1], &s1.vec[2], &s1.vec[3], &s1.vec[4], &s2.vec[0], &s2.vec[1], &s2.vec[2], rhoprime, 0, 1, 2, 3, 4, 5, 6, 7);
        poly_uniform_eta_8x(&s2.vec[3], &s2.vec[4], &s2.vec[5], &t0, &t2, &t3, &t4, &t5, rhoprime, 8, 9, 10, 11, 12, 13, 14, 15);
#elif K == 8 && L == 7
        poly_uniform_eta_8x(&s1.vec[0], &s1.vec[1], &s1.vec[2], &s1.vec[3], &s1.vec[4], &s1.vec[5], &s1.vec[6], &s2.vec[0], rhoprime, 0, 1, 2, 3, 4, 5, 6, 7);
        poly_uniform_eta_8x(&s2.vec[1], &s2.vec[2], &s2.vec[3], &s2.vec[4], &s2.vec[5], &s2.vec[6], &s2.vec[7], &t0, rhoprime, 8, 9, 10, 11, 12, 13, 14, 15);
#else
#error
#endif

        /* Pack secret vectors */
        for (i = 0; i < L; i++)
            polyeta_pack(sk + 3 * SEEDBYTES + i * POLYETA_PACKEDBYTES, &s1.vec[i]);
        for (i = 0; i < K; i++)
            polyeta_pack(sk + 3 * SEEDBYTES + (L + i) * POLYETA_PACKEDBYTES, &s2.vec[i]);

        /* Transform s1 */
        polyvecl_ntt(&s1);

#ifdef DILITHIUM_USE_AES
        aes256ctr_init(&aesctx, rho, 0);
#else
        polyvec_matrix_expand(mat, rho);
#endif

        for (i = 0; i < K; i++)
        {
            /* Expand matrix row */
#ifdef DILITHIUM_USE_AES
            for (unsigned int j = 0; j < L; j++)
            {
                nonce = (i << 8) + j;
                aesctx.n = _mm_loadl_epi64((__m128i *)&nonce);
                poly_uniform_preinit(&row->vec[j], &aesctx);
                poly_nttunpack(&row->vec[j]);
            }
#else
            row = &mat[i];
#endif

            /* Compute inner-product */
            polyvecl_pointwise_acc_montgomery(&t1, row, &s1);
            poly_invntt_tomont(&t1);

            /* Add error polynomial */
            poly_add(&t1, &t1, &s2.vec[i]);

            /* Round t and pack t1, t0 */
            poly_caddq(&t1);
            poly_power2round(&t1, &t0, &t1);
            polyt1_pack(pk + SEEDBYTES + i * POLYT1_PACKEDBYTES, &t1);
            polyt0_pack(sk + 3 * SEEDBYTES + (L + K) * POLYETA_PACKEDBYTES + i * POLYT0_PACKEDBYTES, &t0);
        }

        /* Compute H(rho, t1) and store in secret key */
        shake256(sk + 2 * SEEDBYTES, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
        


        i=0;
        if((memcmp(pk,pk_f,CRYPTO_PUBLICKEYBYTES))!=0||memcmp(sk,sk_f,CRYPTO_SECRETKEYBYTES)!=0)
        {
      
                // for(int k=0;k<256;k++){
                //     printf("%u ",&seedbuf[i]);
                //     if(k%8==7) printf("\n");
                // }
                printf("fail\n");
                return 0;
        }
        else printf("ok\n");
        j++;
    }

    return 0;
}
